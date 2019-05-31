function net = GetNeuronPositions(net)
% set the minimum distance allowed between two cells
limit = 10; % um

% calculate the number of cells in the network
net.MS.N  = round(prod(net.PhysicalDimensions./1000) .* net.CellsPerMillimeterCube); % get the number of MS neurons in the slice
net.FS.N = round(net.MS.N .* (net.FSpercentage/100)); 
N = net.MS.N + net.FS.N;

% get the neuron positions
P = get_positions(net.PhysicalDimensions, N, limit);

% --------------------------------------------------------------
% randomly assign the cells to the MS and FS neuron populations.
x = randperm(N);
net.MS.Position = P(x(1:net.MS.N),:);
net.FS.Position = P(x(net.MS.N+1:end),:);

% randomly assign MSN's to the D1 and D2 class
x = randperm(net.MS.N);
net.MS.D1inds = sort(x(1:round(net.MS.N/2)));
net.MS.D2inds = sort(x((round(net.MS.N/2)+1):net.MS.N));

% -------------------------------------------------------------------------
% get the positions of the cells
function [P] = get_positions(dims, N, limit)
P = rand(N,3);
P(:,1) = P(:,1).*dims(1);
P(:,2) = P(:,2).*dims(2);
P(:,3) = P(:,3).*dims(3);

% move any cell that is closer to its nabours than the limit
P = check_dist(P, dims, limit); 

% -------------------------------------------------------------------------
% set the distance between cells so that it is never less than the limit
function P = check_dist(P, dims, limit)
i = 0;
while i < length(P)
    i = i+1;
    d = get_distance(P, P(i,:));
    inds = find(d<limit);

    selfind = find(inds==i);
    inds = inds([1:selfind-1 selfind+1:end]);
    
    if ~isempty(inds)   
        P(i,:) = rand(1,3) .* dims;
        i = i-1;
    end
end

% -------------------------------------------------------------------------
% calculate the distance between one cell, and all others in the network
function d = get_distance(targetcoords, sourcecoords)
[n,m] = size(targetcoords);
a = targetcoords(:,1) - sourcecoords(1);
b = targetcoords(:,2) - sourcecoords(2);
c = targetcoords(:,3) - sourcecoords(3);
d = sqrt(a.^2 + b.^2 + c.^2);
