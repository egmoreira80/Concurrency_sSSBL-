clear all;
close all;
tic
addpath('common_functions');
addpath('simulation_data');
addpath('concurrency_evaluation');
addpath(genpath('D:\Tools\brainstorm3'));
addpath('D:\Tools\fieldtrip-master');
ft_defaults

%% Human connectome MEG preprocessed file
preproced_data_path = 'E:\Neuroinformatics Collaboratory\Ariosky Areces Gonzalez - Last Simulation\Restin\rmegpreproc\175237_MEG_5-Restin_rmegpreproc.mat';

%% Brainstorm protocol files for MEG/EEG Lead Fields and MEG channel info
protocol_path = 'E:\Neuroinformatics Collaboratory\Ariosky Areces Gonzalez - Last Simulation\BS_protocol';

%% Parameters
method_label  = {'sSSBL' 'sSSBLparcelled' 'sSSBLlaplacian' 'sSSBL2Disotropic' 'sSSBL3Disotropic' 'sSSBLunwrapped' 'sSSBL++' 'eLORETA' 'LCMV'};
% sSSBL with priors of smoothness in 2D cartesian space of samples and band frequencies
% sSSBLparcelled with priors of smoothness in 3D cartesian space of samples, band frequencies and parcells
% sSSBL3Dfield with priors of rotational invariance of 3D source fields 
% sSSBLunwrapped with priors of compensation for surface curvature
% sSSBLplus with all prior information

bands         = [0.1 4; 4.1 8; 8.1 12; 12.1 16; 16.1 32];
band_label    = {'delta' 'theta' 'alpha' 'beta' 'gamma'};

%% Prepare data for MEG/EEG concurrency simulation
sim_data      = struct;

[sim_data]    = meeg_sim_structural(sim_data,preproced_data_path,protocol_path);

[sim_data]    = meeg_sim_functional(sim_data,preproced_data_path);

%% Concurrency evaluation
disp('-->> performing concurrency evaluation');
% outputs: 
% J     - 1xlength(bands) cell array with MEG and EEG source activity vectors concatenated along the second dimension
% stat  - 1xlength(bands) cell array with MEG and EEG source statistic vectors concatenated along the second dimension   
% indms - 1xlength(bands) cell array with MEG and EEG active indexes in concatenated in 1x2 cell array 

%% sSSBL family
ismethod  = 1;

%% sSSBL
disp('-->> Spectral SSBL');
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)

 
[J_sSSBL,stat_sSSBL,indms_sSSBL,data_sSSBL] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_sSSBL,emd_sSSBL] = figure_concurrency(data_sSSBL,sim_data,J_sSSBL,indms_sSSBL,band_label,method_label{1});

%% sSSBLparcelled
disp('-->> Spectral SSBL with parcell smoothness');
isparcel  = 1; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)

 
[J_sSSBLparcelled,stat_sSSBLparcelled,indms_sSSBLparcelled,data_sSSBLparcelled] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_sSSBLparcelled,emd_sSSBLparcelled] = figure_concurrency(data_sSSBLparcelled,sim_data,J_sSSBLparcelled,indms_sSSBLparcelled,band_label,method_label{2});

%% sSSBLlaplacian
disp('-->> Spectral SSBL with laplacian smoothness');
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 1; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
 
[J_sSSBLlaplacian,stat_sSSBLlaplacian,indms_sSSBLlaplacian,data_sSSBLlaplacian] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_sSSBLlaplacian,emd_sSSBLlaplacian] = figure_concurrency(data_sSSBLlaplacian,sim_data,J_sSSBLlaplacian,indms_sSSBLlaplacian,band_label,method_label{3});

%% sSSBL2Disotropy
disp('-->> Spectral SSBL with 2D isotropic rotational invariance');
isparcel  = 0; % 0 (no parcellation smoothness) 1 (parcellation smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 2; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
 
[J_sSSBL2Disotropy,stat_sSSBL2Disotropy,indms_sSSBL2Disotropy,data_sSSBL2Disotropy] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_sSSBL2Disotropy,emd_sSSBL2Disotropy] = figure_concurrency(data_sSSBL2Disotropy,sim_data,J_sSSBL2Disotropy,indms_sSSBL2Disotropy,band_label,method_label{4});

%% sSSBL3Disotropy
disp('-->> Spectral SSBL with 3D isotropic rotational invariance');
isparcel  = 0; % 0 (no parcellation smoothness) 1 (parcellation smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 3; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
 
[J_sSSBL3Disotropy,stat_sSSBL3Disotropy,indms_sSSBL3Disotropy,data_sSSBL3Disotropy] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_sSSBL3Disotropy,emd_sSSBL3Disotropy] = figure_concurrency(data_sSSBL3Disotropy,sim_data,J_sSSBL3Disotropy,indms_sSSBL3Disotropy,band_label,method_label{5});

%% sSSBLunwrapped
disp('-->> Spectral SSBL with unwrapping of curvature');
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)

% giri compensation
iscurv    = 1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_sSSBLgiri,stat_sSSBLgiri,indms_sSSBLgiri,data_sSSBLgiri] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

% sulci compensation
iscurv    = -1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_sSSBLsulc,stat_sSSBLsulc,indms_sSSBLsulc,data_sSSBLsulc] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

% Combining giri and sulci information
J_sSSBLunwrapped     = cell(1,length(bands));
stat_sSSBLunwrapped  = cell(1,length(bands));
indms_sSSBLunwrapped = cell(1,length(bands));
for band = 1:length(bands)
    % --- MEG solutions ---
    data_sSSBLunwrapped{band}{1}   = (data_sSSBLgiri{band}{1}+data_sSSBLsulc{band}{1})/2;
    
    % unwrapped stat
    tmp_stat_giri                  = stat_sSSBLgiri{band}(:,1);
    tmp_stat_sulc                  = stat_sSSBLsulc{band}(:,1);
    tmp_stat                       = (tmp_stat_giri + tmp_stat_sulc)/2;
    stat_sSSBLunwrapped{band}(:,1) = tmp_stat;
    
    % unwrapped indms
    tmp_indms_giri                 = indms_sSSBLgiri{band}{1};
    tmp_indms_sulc                 = indms_sSSBLsulc{band}{1};
    tmp_indms                      = unique([tmp_indms_giri;tmp_indms_sulc]);
    indms_sSSBLunwrapped{band}{1}  = tmp_indms;
    
    % unwrapped J
    tmp_J_giri                     = J_sSSBLgiri{band}(:,1);
    tmp_J_sulc                     = J_sSSBLsulc{band}(:,1);
    tmp_J                          = zeros(length(tmp_stat),1);
    tmp_J(tmp_indms)               = (tmp_J_giri(tmp_indms) + tmp_J_sulc(tmp_indms))/2;
    J_sSSBLunwrapped{band}(:,1)    = tmp_J;
    
    % --- EEG solutions ---
    data_sSSBLunwrapped{band}{2}   = (data_sSSBLgiri{band}{2}+data_sSSBLsulc{band}{2})/2;
    
    % unwrapped stat
    tmp_stat_giri                  = stat_sSSBLgiri{band}(:,2);
    tmp_stat_sulc                  = stat_sSSBLsulc{band}(:,2);
    tmp_stat                       = (tmp_stat_giri + tmp_stat_sulc)/2;
    stat_sSSBLunwrapped{band}(:,2) = tmp_stat;
    
    % unwrapped indms
    tmp_indms_giri                 = indms_sSSBLgiri{band}{2};
    tmp_indms_sulc                 = indms_sSSBLsulc{band}{2};
    tmp_indms                      = unique([tmp_indms_giri;tmp_indms_sulc]);
    indms_sSSBLunwrapped{band}{2}  = tmp_indms;
    
    % unwrapped J
    tmp_J_giri                     = J_sSSBLgiri{band}(:,2);
    tmp_J_sulc                     = J_sSSBLsulc{band}(:,2);
    tmp_J                          = zeros(length(tmp_stat),1);
    tmp_J(tmp_indms)               = (tmp_J_giri(tmp_indms) + tmp_J_sulc(tmp_indms))/2;
    J_sSSBLunwrapped{band}(:,2)    = tmp_J;
end 

[mapsfig_sSSBLunwrapped,emd_sSSBLunwrapped] = figure_concurrency(data_sSSBLunwrapped,sim_data,J_sSSBLunwrapped,indms_sSSBLunwrapped,band_label,method_label{6});

%% sSSBL++
disp('-->> Spectral SSBL with all prior information');
isparcel  = 1; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 1; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)

isfield   = 2; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)

% giri compensation
iscurv    = 1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_sSSBLpp_giri,stat_sSSBLpp_giri,indms_sSSBLpp_giri,data_sSSBLpp_giri] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

% sulci compensation
iscurv    = -1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_sSSBLpp_sulc,stat_sSSBLpp_sulc,indms_sSSBLpp_sulc,data_sSSBLpp_sulc] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

% Combining giri and sulci information
J_sSSBLpp     = cell(1,length(bands));
stat_sSSBLpp  = cell(1,length(bands));
indms_sSSBLpp = cell(1,length(bands));
for band = 1:length(bands)
    % --- MEG solutions ---
    data_sSSBLpp{band}{1}   = (data_sSSBLpp_giri{band}{1}+data_sSSBLpp_sulc{band}{1})/2;
    
    % unwrapped stat
    tmp_stat_giri                  = stat_sSSBLpp_giri{band}(:,1);
    tmp_stat_sulc                  = stat_sSSBLpp_sulc{band}(:,1);
    tmp_stat                       = (tmp_stat_giri + tmp_stat_sulc)/2;
    stat_sSSBLpp{band}(:,1)        = tmp_stat;
    
    % unwrapped indms
    tmp_indms_giri                 = indms_sSSBLpp_giri{band}{1};
    tmp_indms_sulc                 = indms_sSSBLpp_sulc{band}{1};
    tmp_indms                      = unique([tmp_indms_giri;tmp_indms_sulc]);
    indms_sSSBLpp{band}{1}         = tmp_indms;
    
    % unwrapped J
    tmp_J_giri                     = J_sSSBLpp_giri{band}(:,1);
    tmp_J_sulc                     = J_sSSBLpp_sulc{band}(:,1);
    tmp_J                          = zeros(length(tmp_stat),1);
    tmp_J(tmp_indms)               = (tmp_J_giri(tmp_indms) + tmp_J_sulc(tmp_indms))/2;
    J_sSSBLpp{band}(:,1)           = tmp_J;
    
    % --- EEG solutions ---
    data_sSSBLpp{band}{2}   = (data_sSSBLpp_giri{band}{2}+data_sSSBLpp_sulc{band}{2})/2;
    
    % unwrapped stat
    tmp_stat_giri                  = stat_sSSBLpp_giri{band}(:,2);
    tmp_stat_sulc                  = stat_sSSBLpp_sulc{band}(:,2);
    tmp_stat                       = (tmp_stat_giri + tmp_stat_sulc)/2;
    stat_sSSBLpp{band}(:,2)        = tmp_stat;
    
    % unwrapped indms
    tmp_indms_giri                 = indms_sSSBLpp_giri{band}{2};
    tmp_indms_sulc                 = indms_sSSBLpp_sulc{band}{2};
    tmp_indms                      = unique([tmp_indms_giri;tmp_indms_sulc]);
    indms_sSSBLpp{band}{2}         = tmp_indms;
    
    % unwrapped J
    tmp_J_giri                     = J_sSSBLpp_giri{band}(:,2);
    tmp_J_sulc                     = J_sSSBLpp_sulc{band}(:,2);
    tmp_J                          = zeros(length(tmp_stat),1);
    tmp_J(tmp_indms)               = (tmp_J_giri(tmp_indms) + tmp_J_sulc(tmp_indms))/2;
    J_sSSBLpp{band}(:,2)           = tmp_J;
end 

[mapsfig_sSSBLpp,emd_sSSBLpp] = figure_concurrency(data_sSSBLpp,sim_data,J_sSSBLpp,indms_sSSBLpp,band_label,method_label{7});

%% eLORETA family
ismethod  = 2;

%% eLORETA
disp('-->> Spectral eLORETA');
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)

 
[J_eLORETA,stat_eLORETA,indms_eLORETA,data_eLORETA] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_eLORETA,emd_eLORETA] = figure_concurrency(data_eLORETA,sim_data,J_eLORETA,indms_eLORETA,band_label,method_label{8});

%% LCMV family
ismethod  = 3;

%% LCMV
disp('-->> Spectral LCMV');
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)

 
[J_LCMV,stat_LCMV,indms_LCMV,data_LCMV] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

[mapsfig_LCMV,emd_LCMV] = figure_concurrency(data_LCMV,sim_data,J_LCMV,indms_LCMV,band_label,method_label{9});