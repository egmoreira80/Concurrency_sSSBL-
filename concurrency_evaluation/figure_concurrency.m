function [mapsfig,emd] = figure_concurrency(data,sim_data,J,indms,band_label,method_label)
%% Parameters
%% Colormap
load('mycolor.mat');
mapsfig = figure('Color','w','Name',method_label);
surface              = sim_data.structural.surface;
EEG_Channel          = sim_data.structural.EEG_channel.Channel;
MEG_Channel          = sim_data.structural.MEG_channel.Channel;
smoothValue          = 0.66;
SurfSmoothIterations = 10;
surface.Vertices     = tess_smooth(surface.Vertices, smoothValue, SurfSmoothIterations, surface.VertConn, 1);
emd                  = cell(3,6);
emd{2,1}             = 'abs_emd';
emd{3,1}             = 'rel_emd';
emd(1,2:6)           = band_label;

%%
elec_meg             = [];
elec_meg.pos         = zeros(length(MEG_Channel),3);
for ii = 1:length(MEG_Channel)
    elec_meg.lbl{ii}   = MEG_Channel(ii).Name;
    temp               = MEG_Channel(ii).Loc;
    elec_meg.pos(ii,:) = mean(temp,2);
end
elec_meg.label       = elec_meg.lbl;
elec_meg.elecpos     = elec_meg.pos;
elec_meg.unit        = 'mm';
elec_eeg             = [];
elec_eeg.pos = zeros(length(EEG_Channel),3);
for ii = 1:length(EEG_Channel)
    elec_eeg.lbl{ii}   = EEG_Channel(ii).Name;
    temp               = EEG_Channel(ii).Loc;
    elec_eeg.pos(ii,:) = temp;
end
elec_eeg.label       = elec_eeg.lbl;
elec_eeg.elecpos     = elec_eeg.pos;
elec_eeg.unit        = 'mm';

for band = 1:length(band_label)
    %% MEG topography
    subplot(4,5,band)
    temp             = diag(data{band}{1});
    temp             = abs(temp)/max(abs(temp(:)));
    cfg              = [];
    meg              = [];
    cfg.layout       = '4D248_helmet.mat';
    cfg.channel      = 'meg';
    cfg.markers      = '.';
    cfg.markersymbol = '.';
    cfg.colormap     = cmap;
    cfg.markersize   = 3;
    cfg.markercolor  = [1 1 1];
    meg.sens         = elec_meg;
    meg.tra          = elec_meg.pos;
    meg.coilpos      = elec_meg.pos;
    meg.label        = elec_meg.lbl;
    meg.dimord       = 'chan_freq';
    meg.freq         = band;
    meg.powspctrm    = temp;
    ft_topoplotTFR(cfg,meg);
    title(['MEG' ' ' band_label{band} ' ' 'topography']) 
    
    %% MEG source
    subplot(4,5,band+5)
    JMEG                   = J{band}(:,1);
    JMEG                   = JMEG/sum(JMEG);
    mapMEG                 = zeros(length(JMEG),1);
    mapMEG(indms{band}{1}) = JMEG(indms{band}{1});
    patch('Faces',surface.Faces,'Vertices',surface.Vertices,'FaceVertexCData',sim_data.structural.Sulc*0.06+...
        log(1+mapMEG/max(mapMEG)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(-90,90);
    title([method_label ' ' 'MEG' ' ' band_label{band}])
    
    %% EEG topography
    subplot(4,5,band+10)
    temp             = data{band}{2};
    temp             = abs(temp)/max(abs(temp(:)));
    cfg              = [];
    eeg              = [];
    cfg.marker       = '';
    cfg.layout       = 'EEG1020.lay';
    cfg.channel      = 'eeg';
    cfg.markersymbol = '.';
    cfg.colormap     = cmap;
    cfg.markersize   = 3;
    cfg.markercolor  = [1 1 1];
    eeg.elec         = elec_eeg;
    eeg.label        = elec_eeg.lbl;
    eeg.dimord       = 'chan_freq';
    eeg.freq         = band;
    eeg.powspctrm    = temp;
    ft_topoplotTFR(cfg,eeg);
    title(['EEG' ' ' band_label{band} ' ' 'topography'])
    
    %% EEG sources
    subplot(4,5,band+15)
    JEEG                   = J{band}(:,2);
    JEEG                   = JEEG/sum(JEEG);
    mapEEG                 = zeros(length(JEEG),1);
    mapEEG(indms{band}{2}) = JEEG(indms{band}{2});
    patch('Faces',surface.Faces,'Vertices',surface.Vertices,'FaceVertexCData',sim_data.structural.Sulc*0.06+...
        log(1+mapEEG/max(mapEEG)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    set(gca,'xcolor','w','ycolor','w','zcolor','w'),
    view(-90,90);
    title([method_label ' ' 'EEG' ' ' band_label{band}])
    
    %% Quality measures
    [emd_val]        = estimating_EMD(mapMEG,mapEEG,surface.Vertices,surface.Faces);
    emd{2,band+1}    = emd_val;
end

JMEGdelta                   = J{1}(:,1);
JMEGdelta                   = JMEGdelta/sum(JMEGdelta);
mapMEGdelta                 = zeros(length(JMEGdelta),1);
mapMEGdelta(indms{1}{1})    = JMEGdelta(indms{1}{1});

JMEGalpha                   = J{3}(:,1);
JMEGalpha                   = JMEGalpha/sum(JMEGalpha);
mapMEGalpha                 = zeros(length(JMEGalpha),1);
mapMEGalpha(indms{3}{1})    = JMEGalpha(indms{3}{1});

[emd_val]                   = estimating_EMD(mapMEGdelta,mapMEGalpha,surface.Vertices,surface.Faces);
for band = 1:length(band_label)
    emd{3,band+1}    = emd{2,band+1}/emd_val;
end

end
