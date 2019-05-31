function [h] = raster_plot(events,times,varargin)

% RASTER_PLOT raster plot
%   RASTER_PLOT(E,T) plots the spike events in E at corresponding times T as a raster plot, with one point per spike, one row
%   per event index (i.e. either per neuron or per sweep)
%
%   RASTER_PLOT(E,T,FLAG) where FLAG is:
%       'r': randomises the order in which the events are plotted - this is useful for removing any potentially arbitrary structure
%       imposed by the order of event indexing (e.g. in the BG models, the events are ordered by channel)
%
%       's': puts the rasterplot as the top window of a 2x1 subplot [put 'rs' to get both]
%
%   RASTER_PLOT(E,T,FLAG,STRING) adds the STRING as the title of the plot (put FLAG=[] to omit)
%
%   Returns the handle to the figure window
%
%   Mark Humphries 9/10/2009


new_events = events;
if nargin >= 3 & findstr(varargin{1},'r')
    % new_times = [];
    event_idxs = unique(events);               % array of indices
    rand_seq = randperm(length(event_idxs));
    map = event_idxs(rand_seq);                % array of indices to re-map to
   
    for loop=1:length(map)
        new_events(events==event_idxs(loop)) = map(loop);   % replace     
    end
end

h = figure 
if nargin >= 3 & findstr(varargin{1},'s')
    subplot(211)
end
plot(times,new_events,'k.')
min(new_events);
%axis([min(new_times) max(new_times) min(new_events) max(new_events)]);

if nargin==4
    title(varargin{2});
end
ylabel('event')
xlabel('time');



