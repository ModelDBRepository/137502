%% testing striatal model construction using "fully-wired"
% approximation....

clear all

%% network parameters
SIMPARAMS.net.PhysicalDimensions = [300 300 300];   % slice dimensions
SIMPARAMS.net.CellsPerMillimeterCube = 84900;       % number of MSNs per mm^3 (Oorschot: 84,900 MSNs per mm3 of rat)

SIMPARAMS.net.FSpercentage = 3;                     % percentage of FS neurons (N_MS * FSpercentage/100)
% "fully-wired" model, target numbers taken from striatal anatomical network model paper
SIMPARAMS.net.ConnectMethod = 'random';     % either 'physical' or 'random'. 
SIMPARAMS.net.random.targetN_msms =  422;   % target number of MSNS INPUT TO 1 MSN connections for this % FSI..
SIMPARAMS.net.random.targetN_fsms =  7.88;   % target number of FSIS INPUT TO 1 MSN connections for this % FSI..
SIMPARAMS.net.random.targetN_fsfs =  3.29;     % target number of FSIS INPUT TO 1 FSI connections for this % FSI..
SIMPARAMS.net.random.targetN_fsgap =  2.04;  % target number of gap junctions per FSI for this % of FSIs

% SIMPARAMS.net.FSpercentage = 5;                     % percentage of FS neurons (N_MS * FSpercentage/100)
% % "fully-wired" model, target numbers taken from striatal anatomical network model paper
% SIMPARAMS.net.ConnectMethod = 'random';     % either 'physical' or 'random'. 
% SIMPARAMS.net.random.targetN_msms =  422;   % target number of MSNS INPUT TO 1 MSN connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsms =  13.3;   % target number of FSIS INPUT TO 1 MSN connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsfs =  5.95;     % target number of FSIS INPUT TO 1 FSI connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsgap =  3.72;  % target number of gap junctions per FSI for this % of FSIs


SIMPARAMS.net.ConnectMethod = 'physical'; 

% set the physiological parameters for the intrastriatal pathways
% required to fill in parameters that the striatum building function
% calls...
SIMPARAMS.physiology.MS_MS.maxdelay = 2; % msec
SIMPARAMS.physiology.MS_MS.baseweight = 4.36;
SIMPARAMS.physiology.MS_MS.baseweightSD = 0;
SIMPARAMS.physiology.FS_MS.maxdelay = 2; % msec
SIMPARAMS.physiology.FS_MS.baseweight = (4.36 * 5);
SIMPARAMS.physiology.FS_MS.baseweightSD = 0;
SIMPARAMS.physiology.FS_FS.maxdelay = 2; % msec
SIMPARAMS.physiology.FS_FS.baseweight = (4.36 * 5);
SIMPARAMS.physiology.FS_FS.baseweightSD = 0;
SIMPARAMS.physiology.FS_gap.baseweight = (150 / 5);
SIMPARAMS.physiology.FS_gap.baseweightSD = 0;

%%% set some simulation parameters required by construction function...
SIMPARAMS.sim.RANDSEED = int32(12345);
SIMPARAMS.sim.logfname = 'logfile.log';

SIMPARAMS.sim.tstart = 0;           % msec : double
SIMPARAMS.sim.tfinal = 100;         % msec : double
SIMPARAMS.sim.dt = 0.1;            % msec : double

%% run code 
% get positions, and number of MSNs and FSIs
SIMPARAMS.net = GetNeuronPositions(SIMPARAMS.net); 

% wire them up...
net = BuildStriatumNetwork(SIMPARAMS.net,SIMPARAMS.physiology,SIMPARAMS.sim);

% shift indices
Cmsms = net.Cmsms +1; % was zero-indexed
Cmsms_b = net.Cmsms_b+1; % was zero-indexed
Cfsms = net.Cfsms +1; % was zero-indexed
Cfsms_b = net.Cfsms_b+1; % was zero-indexed
Cfsfs = net.Cfsfs +1; % was zero-indexed
Cfsfs_b = net.Cfsfs_b+1; % was zero-indexed
Pgapfs = net.Pgapfs +1; % was zero-indexed

% storage
nMSN_MSNtgts = []; MSN_MSNtgt_d = {}; 
nMSN_MSNsrcs = []; MSN_MSNsrc_d = {}; 
nFSI_MSNsrcs = []; FSI_MSNsrc_d = {}; 
nFSI_MSNtgts = []; FSI_MSNtgt_d = {}; 
nFSI_FSItgts = []; FSI_FSItgt_d = {}; 
nFSI_FSIsrcs = []; FSI_FSIsrc_d = {}; 
nFSIgap = []; FSIgap_d = {}; 

%% MSN connections
for idx = 1:SIMPARAMS.net.MS.N
    thispos = net.MS.Position(idx,:);
    %%% output to other MSNs
    tgts = Cmsms(Cmsms_b(idx):Cmsms_b(idx+1)-1);
    nMSN_MSNtgts = [nMSN_MSNtgts; numel(tgts)];     % build array: number of MSN tgts
    tgtpos = net.MS.Position(tgts,:);   % position of all targets
    % array of all distances to tgts
    dists = sqrt(sum((tgtpos - repmat(thispos,numel(tgts),1))'.^2))';
    MSN_MSNtgt_d = [MSN_MSNtgt_d; dists ]; % distances to all targets

    %%% input
    % inputs from other MSN
    MSNbs = find(Cmsms == idx); % idx of srcs
    MSNsrcs = zeros(numel(MSNbs),1);
    for k = 1:numel(MSNbs)
        % the source cell idx corresponds to the idx of the bound array
        % entry lower than the idx in Cmsms
        MSNsrcs(k) = find(Cmsms_b <= MSNbs(k),1,'last'); 
    end
    nMSN_MSNsrcs = [nMSN_MSNsrcs; numel(MSNsrcs)];  % number of source MSNs
    srcpos = net.MS.Position(MSNsrcs,:);   % position of all sources
    dists = sqrt(sum((srcpos - repmat(thispos,numel(MSNsrcs),1))'.^2))'; % distances to all targets
    MSN_MSNsrc_d = [MSN_MSNsrc_d; dists];

    % inputs from FSIs
    FSIbs = find(Cfsms == idx); % idx of srcs
    FSIsrcs = zeros(numel(FSIbs),1);
    for k = 1:numel(FSIbs)
        % the source cell idx corresponds to the idx of the bound array
        % entry lower than the idx in Cmsms
        FSIsrcs(k) = find(Cfsms_b <= FSIbs(k),1,'last'); 
    end
    nFSI_MSNsrcs = [nFSI_MSNsrcs; numel(FSIsrcs)];  % build array: number of source FSIs
    srcpos = net.FS.Position(FSIsrcs,:);   % position of all sources
    dists = sqrt(sum((srcpos - repmat(thispos,numel(FSIsrcs),1))'.^2))'; % distances to all targets
    FSI_MSNsrc_d = [FSI_MSNsrc_d; dists]; % cell array of individual distributions
end

%% FSI connection stats
for idx = 1:SIMPARAMS.net.FS.N
    %%% output to MSNs
    tgts = Cfsms(Cfsms_b(idx):Cfsms_b(idx+1)-1);
    tgtpos = net.MS.Position(tgts,:);   % position of all targets
    % array of all distances to tgts
    dists = sqrt(sum((tgtpos - repmat(thispos,numel(tgts),1))'.^2))'; % distances to all targets
    FSI_MSNtgt_d = [FSI_MSNtgt_d; dists];
    nFSI_MSNtgts = [nFSI_MSNtgts; numel(tgts)];     % build array: number of MSN tgts

    % to other FSIs
    tgts = Cfsfs(Cfsfs_b(idx):Cfsfs_b(idx+1)-1);
    tgtpos = net.FS.Position(tgts,:);   % position of all targets
    % array of all distances to tgts - note that this is stored so
    dists = sqrt(sum((tgtpos - repmat(thispos,numel(tgts),1))'.^2))'; % distances to all targets
    FSI_FSItgt_d = [FSI_FSItgt_d; dists];
    nFSI_FSItgts = [nFSI_FSItgts; numel(tgts)];     % build array: number of MSN tgts

    %%% input from other FSIs (axo-dendritic)
    FSIbs = find(Cfsfs == idx); % idx of srcs
    FSIsrcs = zeros(numel(FSIbs),1);
    for k = 1:numel(FSIbs)
        % the source cell idx corresponds to the idx of the bound array
        % entry lower than the idx in Cmsms
        FSIsrcs(k) = find(Cfsfs_b <= FSIbs(k),1,'last'); 
    end
    srcpos = net.FS.Position(FSIsrcs,:);   % position of all sources
    dists = sqrt(sum((srcpos - repmat(thispos,numel(FSIsrcs),1))'.^2))'; % distances to all targets
    FSI_FSIsrc_d = [FSI_FSIsrc_d; dists]; 
    nFSI_FSIsrcs = [nFSI_FSIsrcs; numel(FSIsrcs)];  % build array

    %%% FS-FS gap junctions....
    % find all gap junctions involving the current FSI
    c1 = find(Pgapfs(:,1) == idx); c2 = find(Pgapfs(:,2) == idx);
    gapvec = [c1; c2];
    gaplist = [fliplr(Pgapfs(c2,:)); Pgapfs(c1,:)]; % create list of all gap junctions
    gappos = net.FS.Position(gaplist(:,2),:);   % position of all gap junction coupled neurons
    dists = sqrt(sum((gappos - repmat(thispos,numel(c1)+numel(c2),1))'.^2))'; % distances to all targets
    FSIgap_d = [FSIgap_d; dists];
    nFSIgap = [nFSIgap; numel(gapvec)]; % vector of number of gap junctions

end

