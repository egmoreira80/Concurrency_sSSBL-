function [sim_data] = meeg_sim_functional(sim_data,preproced_data_path)
%% Loading data
preproced_data = load(preproced_data_path);
%%
disp('-->> Creating cross-spectrum');
Ntpoints    = size(preproced_data.data.trial{1,1},2);
Nsegments   = length(preproced_data.data.trial);
Nsens       = length(preproced_data.data.label);
Nw          = 3;
Data_FC     = zeros(Nsens,Ntpoints,2*Nw,Nsegments);
Fs          = preproced_data.data.fsample;
deltaf      = Fs/Ntpoints;
F           = 0:deltaf:((Ntpoints-1)*deltaf);
Nfreq       = length(F);
Svv         = zeros(Nsens,Nsens,Nfreq);
e           = dpss(Ntpoints,Nw);
e           = reshape(e,[1,Ntpoints,2*Nw]);
for seg = 1:Nsegments
    disp(strcat('-->> processing time segment: ', num2str(seg)));
    tmp                   = preproced_data.data.trial{1,seg};
    tmp                   = repmat(tmp,[1,1,2*Nw]).*repmat(e,[Nsens,1,1]);
    tmp_FC                = fft(tmp,[],2);
    Data_FC(:,:,:,seg)    = tmp_FC;
    for freq = 1:Nfreq
        Svv(:,:,freq)     = Svv(:,:,freq) + squeeze(tmp_FC(:,freq,:))*squeeze(tmp_FC(:,freq,:))';
    end
end
Svv = Svv/Nsegments;

%% Saving data
sim_data.functional.F            = F;
sim_data.functional.Svv          = Svv;
sim_data.functional.Data_FC      = Data_FC;
sim_data.functional.deltaf       = deltaf;
sim_data.functional.Ntpoints     = Ntpoints;
sim_data.functional.Nsegments    = Nsegments;
end