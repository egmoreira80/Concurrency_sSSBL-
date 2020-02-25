function [L_MEG,MEG_channel] = remove_leadfield_channel(MEG_channel,unproced,preproced)

Channels = MEG_channel.Channel;
labels = preproced.data.label;
L_MEG = unproced.Gain;

from = 1;
limit = length(Channels);
while(from <= limit)
    pos = find(strcmpi(Channels(from).Name, labels), 1);
    if (isempty(pos))
        Channels(from)=[];
        L_MEG(from,:)=[];
        limit = limit - 1;
    else
        from = from + 1;
    end
end

MEG_channel.Channel = Channels;
end

