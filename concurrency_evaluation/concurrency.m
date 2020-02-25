function [J,stat,indms,data] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv)
%% Parameters
F            = sim_data.functional.F;
Svv          = sim_data.functional.Svv;
Nsegments    = sim_data.functional.Nsegments;

%% parcel/field options
if isparcel == 0
    if (isfield == 1) || (isfield == 2)
        parcellation = sim_data.structural.parcellation_none;
    elseif isfield == 3
        parcellation = sim_data.structural.parcellation_none3D;
    end
elseif isparcel == 1
    if (isfield == 1) || (isfield == 2)
        parcellation = sim_data.structural.parcellation;
    elseif isfield == 3
        parcellation = sim_data.structural.parcellation3D;
    end
end

%% neigh/field options
if isneigh == 0
    if isfield == 1
        W  = speye(length(sim_data.structural.D));
    elseif isfield == 2
        W  = sim_data.structural.Ninv;
    elseif isfield == 3
        W  = speye(length(sim_data.structural.D3D));
    end
elseif isneigh == 1
    if isfield == 1
        W  = sim_data.structural.Dinv;
    elseif isfield == 2
        W  = sim_data.structural.DNinv;
    elseif isfield == 3
        W  = sim_data.structural.D3Dinv;
    end
end

%% curv/field options
if iscurv == 0
    if isfield == 1
        L_MEG = sim_data.structural.L_MEG;
        L_EEG = sim_data.structural.L_EEG;
    elseif (isfield == 2) || (isfield == 3)
        L_MEG = sim_data.structural.L_MEG3D;
        L_EEG = sim_data.structural.L_EEG3D;
    end
elseif iscurv == 1
    if isfield == 1
        L_MEG = sim_data.structural.L_MEGgiri;
        L_EEG = sim_data.structural.L_EEGgiri;
    elseif (isfield == 2) || (isfield == 3)
        L_MEG = sim_data.structural.L_MEG3Dgiri;
        L_EEG = sim_data.structural.L_EEG3Dgiri;
    end
elseif iscurv == -1
    if isfield == 1
        L_MEG = sim_data.structural.L_MEGsulc;
        L_EEG = sim_data.structural.L_EEGsulc;
    elseif (isfield == 2) || (isfield == 3)
        L_MEG = sim_data.structural.L_MEG3Dsulc;
        L_EEG = sim_data.structural.L_EEG3Dsulc;
    end
end

%%
for band = 1:length(bands)
    F1          = bands(band,1); %frequency band lower limit
    F2          = bands(band,2); %frequency band upper limit
    [val,idf1]  = min(abs(F-F1));
    [val,idf2]  = min(abs(F-F2));
    SvvMEG      = mean(Svv(:,:,idf1:idf2),3);
    Nsamples    = Nsegments*(idf2-idf1+1);
    
    %% MEG Inverse solution filter
    [T_MEG,J_MEG,statMEG,indmsMEG] = inverse_solver(L_MEG,parcellation,W,SvvMEG,Nsamples,ismethod,isfield);
    
    %% Creating EEG Cross-spectra
    SjjMEG      = T_MEG*SvvMEG*T_MEG';
    SjjMEGth    = zeros(size(T_MEG,1),size(T_MEG,1));   
    if isfield == 1        
        SjjMEGth(indmsMEG,indmsMEG) = SjjMEG(indmsMEG,indmsMEG);
    elseif (isfield == 2) || (isfield == 3)
        for ii = 1:length(indmsMEG)
            q0                                      = 3*indmsMEG(ii);
            SjjMEGth([q0-2;q0-1;q0],[q0-2;q0-1;q0]) = SjjMEG([q0-2;q0-1;q0],[q0-2;q0-1;q0]);
        end
    end
    SvvEEG      = L_EEG*SjjMEGth*L_EEG';
            
    %% EEG inverse solution filter
    [T_EEG,J_EEG,statEEG,indmsEEG] = inverse_solver(L_EEG,parcellation,W,SvvEEG,Nsamples,ismethod,isfield);
    
    %% Saving solutions
    J{band}        = [J_MEG J_EEG];
    stat{band}     = [statMEG statEEG];
    indms{band}    = {indmsMEG indmsEEG}; 
    data{band}     = {SvvMEG SvvEEG};
end