function [sim_data] = meeg_sim_structural(sim_data,preproced_data_path,protocol_path)
%% Loading files
preproced_data = load(preproced_data_path);

leadfield_MEG_file = fullfile(protocol_path,'data','175237_MEG_4D_Neuroimaging','@raw5-Restin_c_rfDC','headmodel_surf_os_meg');
leadfield_MEG = load(leadfield_MEG_file);

leadfield_EEG_file = fullfile(protocol_path,'data','175237_EEG_10_20','@intra','headmodel_surf_openmeeg');
leadfield_EEG = load(leadfield_EEG_file);

MEG_channel_file = fullfile(protocol_path,'data','175237_MEG_4D_Neuroimaging','@raw5-Restin_c_rfDC','channel_4d_acc1.mat');
MEG_channel = load(MEG_channel_file);

EEG_channel_file = fullfile(protocol_path,'data','175237_EEG_10_20','@intra','channel_10-20_19.mat');
EEG_channel = load(EEG_channel_file);

surface_file = fullfile(protocol_path,'anat','175237_MEG_4D_Neuroimaging','tess_cortex_concat_6000V.mat');
surface = load(surface_file);

%% Surface curvature compensator
disp('-->> Creating curvature compensator');
aSulc                 = 5; % baseline of sulci curvature factor 
aGiri                 = 5; % baseline of giri curvature factor 
bSulc                 = 3; % scale of sulci curvature factor 
bGiri                 = 3; % scale of giri curvature factor 

Curv                  = surface.Curvature;
Sulc                  = surface.SulciMap;
Curv                  = abs(Curv);
CurvSulc              = zeros(length(Curv),1);
CurvGiri              = zeros(length(Curv),1);
CurvSulc(Sulc == 1)   = aSulc + bSulc.*Curv(Sulc == 1);
CurvSulc(Sulc == 0)   = 1;
CurvGiri(Sulc == 0)   = aGiri + bGiri.*Curv(Sulc == 0);
CurvGiri(Sulc == 1)   = 1;

Sulc3D                = zeros(1,3*length(Sulc));
CurvSulc3D            = zeros(1,3*length(Curv));
CurvGiri3D            = zeros(1,3*length(Curv));

node3 = 1;
for node = 1:length(Curv)
    CurvSulc3D([node3 node3+1 node3+2]) = repmat(CurvSulc(node),1,3);
    CurvGiri3D([node3 node3+1 node3+2]) = repmat(CurvGiri(node),1,3);
    Sulc3D([node3 node3+1 node3+2])     = repmat(Sulc(node),1,3);
    node3                               = node3 + 3;
end

% Compensated Lead Fields
disp('-->> Creating field information');
[L_MEG3D,MEG_channel] = remove_leadfield_channel(MEG_channel,leadfield_MEG,preproced_data);
L_EEG3D               = leadfield_EEG.Gain;

GridOrient            = leadfield_MEG.GridOrient;
GridAtlas             = leadfield_MEG.GridAtlas;

L_MEG                 = bst_gain_orient(L_MEG3D, GridOrient,GridAtlas);
L_EEG                 = bst_gain_orient(L_EEG3D, GridOrient,GridAtlas);

clearvars leadfield_MEG leadfield_EEG

L_MEG3Dsulc           = L_MEG3D.*repmat(CurvSulc3D,size(L_MEG3D,1),1);
L_EEG3Dsulc           = L_EEG3D.*repmat(CurvSulc3D,size(L_EEG3D,1),1);
L_MEG3Dgiri           = L_MEG3D.*repmat(CurvGiri3D,size(L_MEG3D,1),1);
L_EEG3Dgiri           = L_EEG3D.*repmat(CurvGiri3D,size(L_EEG3D,1),1);

L_MEGsulc             = bst_gain_orient(L_MEG3Dsulc,GridOrient,GridAtlas);
L_EEGsulc             = bst_gain_orient(L_EEG3Dsulc,GridOrient,GridAtlas);
L_MEGgiri             = bst_gain_orient(L_MEG3Dgiri,GridOrient,GridAtlas);
L_EEGgiri             = bst_gain_orient(L_EEG3Dgiri,GridOrient,GridAtlas);

%% Parcells 
disp('-->> Creating parcel smoother');
parcellation_none   = cell(length(L_MEG),1);
for area = 1:length(L_MEG)
   parcellation_none{area}    = area;
end

parcellation_none3D = cell(length(L_MEG),1);
for area = 1:length(L_MEG)
    q0                        = 3*(area-1);
    parcellation_none3D{area} = [q0+1;q0+2;q0+3];
end

Atlas = surface.Atlas(3).Scouts;

parcellation        = cell(length(Atlas),1);
for area = 1:length(Atlas)
    parcellation{area}        = Atlas(area).Vertices;
end

parcellation3D      = cell(length(Atlas),1);
for area = 1:length(Atlas)
    for node = 1:length(Atlas(area).Vertices)
        q0                    = 3*(Atlas(area).Vertices(node)-1);
        tmp_parcellation3D    = [q0+1;q0+2;q0+3];
        parcellation3D {area} = cat(1,parcellation3D {area},tmp_parcellation3D );
    end
end

%% Laplacian and Normals
disp('-->> Creating Laplacian&Normals');
Faces   = surface.Faces; 
[D,D3D] = graph_laplacian(Faces);
Dinv    = inv(D);
Dinv    = (Dinv + Dinv)/2;
D3Dinv  = inv(D3D);
D3Dinv  = (D3Dinv + D3Dinv)/2;
Ninv    = blk_diag(GridOrient',1);
N       = Ninv';
DN      = D*N;
DNinv   = Ninv*Dinv;

%% Saving data
sim_data.structural.L_MEG3D             = L_MEG3D;
sim_data.structural.L_MEG3Dgiri         = L_MEG3Dgiri;
sim_data.structural.L_MEG3Dsulc         = L_MEG3Dsulc;
sim_data.structural.L_EEG3D             = L_EEG3D;
sim_data.structural.L_EEG3Dgiri         = L_EEG3Dgiri;
sim_data.structural.L_EEG3Dsulc         = L_EEG3Dsulc;
sim_data.structural.L_MEG               = L_MEG;
sim_data.structural.L_MEGgiri           = L_MEGgiri;
sim_data.structural.L_MEGsulc           = L_MEGsulc;
sim_data.structural.L_EEG               = L_EEG;
sim_data.structural.L_EEGgiri           = L_EEGgiri;
sim_data.structural.L_EEGsulc           = L_EEGsulc;
sim_data.structural.surface             = surface;
sim_data.structural.MEG_channel         = MEG_channel;
sim_data.structural.EEG_channel         = EEG_channel;
sim_data.structural.parcellation_none   = parcellation_none;
sim_data.structural.parcellation_none3D = parcellation_none3D;
sim_data.structural.parcellation        = parcellation;
sim_data.structural.parcellation3D      = parcellation3D;
sim_data.structural.GridOrient          = GridOrient;
sim_data.structural.GridAtlas           = GridAtlas;
sim_data.structural.Sulc                = Sulc;
sim_data.structural.CurvSulc            = CurvSulc;
sim_data.structural.CurvGiri            = CurvGiri;
sim_data.structural.DN                  = DN;
sim_data.structural.DNinv               = DNinv;
sim_data.structural.D                   = D;
sim_data.structural.Dinv                = Dinv;
sim_data.structural.D3D                 = D3D;
sim_data.structural.D3Dinv              = D3Dinv;
sim_data.structural.N                   = N;
sim_data.structural.Ninv                = Ninv;
