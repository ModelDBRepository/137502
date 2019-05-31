% -------------------------------------------------------------------------
% All default parameters for the striatal network model
% -------------------------------------------------------------------------
function SIMPARAMS = StriatumNetworkParameters(varargin)

% pass in argument with full path and filename to load an exisiting
% network.

% set the path to the network building functions
addpath ./CreateNetwork/
blnLoad = 0;
if nargin >= 1 
    blnLoad = 1; 
    network_fname = varargin{1};
end

% network parameters
% SIMPARAMS.net.PhysicalDimensions = [1000 1000 1000];   % slice dimensions for making spatial inputs within series of shells...
SIMPARAMS.net.PhysicalDimensions = [500 500 500];   % slice dimensions for making spatial inputs within single sphere...
% SIMPARAMS.net.PhysicalDimensions = [250 250 250];   % slice dimensions
% SIMPARAMS.net.PhysicalDimensions = [181 181 181];   % slice dimensions to get ~500 MSNs...
% SIMPARAMS.net.PhysicalDimensions = [106 106 106];   % slice dimensions to get ~100 MSNs...
SIMPARAMS.net.CellsPerMillimeterCube = 84900;       % number of MSNs per mm^3
SIMPARAMS.net.FSpercentage = 1;                     % percentage of FS neurons (N_MS * FSpercentage/100)
SIMPARAMS.net.ConnectMethod = 'physical';

% Figures from PLoS Comp Biol paper v1...
% SIMPARAMS.net.ConnectMethod = 'random';     % either 'physical' or 'random'. 
% SIMPARAMS.net.random.targetN_msms =  422;   % target number of MSNS INPUT TO 1 MSN connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsms =  7.88;   % target number of FSIS INPUT TO 1 MSN connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsfs =  3.29;     % target number of FSIS INPUT TO 1 FSI connections for this % FSI..
% SIMPARAMS.net.random.targetN_fsgap =  2.04;  % target number of gap junctions per FSI for this % of FSIs

% % Figures to get Ponzi-like network
% SIMPARAMS.net.ConnectMethod = 'random';     % either 'physical' or 'random'. 
% %SIMPARAMS.net.random.targetN_msms =  265;   % target number of MSNS INPUT TO 1 MSN connections for 250-on-side PW network
% %SIMPARAMS.net.random.targetN_msms =  100;   % target number of MSNS INPUT TO 1 MSN connections for 181-on-side PW network
% SIMPARAMS.net.random.targetN_msms =  20;   % target number of MSNS INPUT TO 1 MSN connections for 106-on-side PW network (p=0.2)
% SIMPARAMS.net.random.targetN_fsms =  1;   % irrelevant as only MSNs are connected
% SIMPARAMS.net.random.targetN_fsfs =  1;     % irrelevant as only MSNs are connected
% SIMPARAMS.net.random.targetN_fsgap =  1;  % irrelevant as only MSNs are connected

% Use the PhysicalDimensions to get MS.N, FS.N, D1inds, D2inds, and the
% neuron positions MS_P, and FS_P. *** NOTE: was CellsPerMillimeterCube
% based on just MSNs, or all cell types? ***
[SIMPARAMS.net] = GetNeuronPositions(SIMPARAMS.net); 

% Assign neurons to the two default channels
SIMPARAMS.net.CHAN1_MS = int32(0:(round(SIMPARAMS.net.MS.N / 2)-1))';
SIMPARAMS.net.CHAN1_FS = int32(0:(round(SIMPARAMS.net.FS.N / 2)-1))';
SIMPARAMS.net.CHAN2_MS = int32(round(SIMPARAMS.net.MS.N / 2):(SIMPARAMS.net.MS.N-1))';
SIMPARAMS.net.CHAN2_FS = int32(round(SIMPARAMS.net.FS.N / 2):(SIMPARAMS.net.FS.N-1))';

% -------------------------------------------------------------------------
% neuron parameters
SIMPARAMS.physiology.MSparams = getMSparams(SIMPARAMS);   % get the parameters for the DA MSN model
SIMPARAMS.physiology.FSparams = getFSparams(SIMPARAMS);   % get the parameters for the DA FSI model 

SIMPARAMS.physiology.Eglu = 0;             % msec : double
SIMPARAMS.physiology.Egaba = -60;          % msec : double
SIMPARAMS.physiology.ts_glu_AMPA = 6;      % msec : double
SIMPARAMS.physiology.ts_glu_NMDA = 160;    % msec : double
SIMPARAMS.physiology.ts_gaba = 4;          % msec : double - use 50ms for Ponzi-like network
SIMPARAMS.physiology.tau_fsgap = 11;       % msec : double

SIMPARAMS.physiology.glu_ratio = 0.5;  % based on DA paper, which was based on Moyer et al  
SIMPARAMS.physiology.DA = 0.0; 

% set the physiological parameters for the intrastriatal pathways
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

% -------------------------------------------------------------------------
% Simulation parameters
SIMPARAMS.sim.RANDSEED = int32(12345);
SIMPARAMS.sim.logfname = 'logfile.log';

SIMPARAMS.sim.tstart = 0;           % msec : double
SIMPARAMS.sim.tfinal = 100;         % msec : double
SIMPARAMS.sim.dt = 0.1;            % msec : double

% load a network from previous file, or build a new one 
if blnLoad
    SIMPARAMS.net = LoadStriatumNetwork(network_fname);
else
    SIMPARAMS.net = BuildStriatumNetwork(SIMPARAMS.net, SIMPARAMS.physiology, SIMPARAMS.sim);
end

SIMPARAMS.sim.RecordChan_MS = int32(0:49)';

SIMPARAMS.sim.MSspikebuffer = int32(2000000); % int32
SIMPARAMS.sim.FSspikebuffer = int32(2000000); % int32

SIMPARAMS.sim.initVms = ones(SIMPARAMS.net.MS.N, 1) .*  SIMPARAMS.physiology.MSparams(:,2);        % mV
SIMPARAMS.sim.initUms = zeros(SIMPARAMS.net.MS.N, 1); 
SIMPARAMS.sim.initVfs = ones(SIMPARAMS.net.FS.N, 1) .*  SIMPARAMS.physiology.FSparams(:,2);        % mV
SIMPARAMS.sim.initUfs = zeros(SIMPARAMS.net.FS.N, 1); 
SIMPARAMS.sim.initVgapfs = ones(length(SIMPARAMS.net.Pgapfs),1) .* SIMPARAMS.physiology.MSparams(1,2); % mV

SIMPARAMS.sim.SpikeEventQue_MS = zeros(SIMPARAMS.net.MS.N, (SIMPARAMS.physiology.MS_MS.maxdelay / SIMPARAMS.sim.dt)+1);
SIMPARAMS.sim.SpikeEventQue_FS = zeros(SIMPARAMS.net.FS.N, (SIMPARAMS.physiology.FS_FS.maxdelay / SIMPARAMS.sim.dt)+1);

SIMPARAMS.sim.Iinj_MS = zeros(SIMPARAMS.net.MS.N,1); 
SIMPARAMS.sim.Iinj_FS = zeros(SIMPARAMS.net.FS.N,1);

% -------------------------------------------------------------------------
% unstructured cortical input
SIMPARAMS.input.CTX.r_MSSEG = ones(SIMPARAMS.net.MS.N,1) .* 1.9; 
SIMPARAMS.input.CTX.N_MSSEG = int32(ones(SIMPARAMS.net.MS.N,1) .* 250);
SIMPARAMS.input.CTX.alpha_MSSEG = ones(SIMPARAMS.net.MS.N,1) .* 0.0;
SIMPARAMS.input.CTX.r_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* 1.9;
SIMPARAMS.input.CTX.N_FSSEG = int32(ones(SIMPARAMS.net.FS.N,1) .* 250);
SIMPARAMS.input.CTX.alpha_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* 0.0;

% The selection experiment parameters
SIMPARAMS.input.Selection.Pt =  int32([-1] ./ SIMPARAMS.sim.dt); % time of change in msec, converted into iteration
SIMPARAMS.input.Selection.Pch = int32([1]'); % channel
SIMPARAMS.input.Selection.Phz = SIMPARAMS.input.CTX.r_MSSEG(1); % rate of spkgen

% parameters for the Nisenbaum cortical pulse experiments
SIMPARAMS.input.PULSE.P = zeros(SIMPARAMS.sim.tfinal.*100,1);
SIMPARAMS.input.PULSE.r_ctx = 0; 
SIMPARAMS.input.PULSE.Nctx_ms = 0; 
SIMPARAMS.input.PULSE.Nctx_fs = 0; 
SIMPARAMS.input.PULSE.ts_spks = 50; % msec

SIMPARAMS.input.PULSE.pulsetimes = [250 750 1250 1750 2250 2750 3250 3750 4250 4750];
SIMPARAMS.input.PULSE.ISI = 0;
SIMPARAMS.input.PULSE.firstpulse = 0;
SIMPARAMS.input.PULSE.secondpulse = SIMPARAMS.input.PULSE.firstpulse + SIMPARAMS.input.PULSE.ISI;
% -------------------------------------------------------------------------
% parameters for the MS neurons - from Humphries et al (2009) NN paper
% -------------------------------------------------------------------------
function MSparams = getMSparams(SIMPARAMS)
MSparams = ones(SIMPARAMS.net.MS.N,14);

MSparams(:,1) = 50;       % C
MSparams(:,2) = -80;      % vr... -66 is approx Ponzi cells...
MSparams(:,3) = -33.8;    % vt 
MSparams(:,4) = 0.05;     % a
MSparams(:,5) = -20;      % b
MSparams(:,6) = -55;      % c
MSparams(:,7) = 377;      % d
MSparams(:,8) = 40;       % vp
MSparams(:,9) = 1.14;     % k
MSparams(:,10) = -68.4;   % ED_ms
MSparams(:,11) = 0.03;    % alpha_ms
MSparams(SIMPARAMS.net.MS.D1inds,11) = MSparams(SIMPARAMS.net.MS.D1inds,11) .* 0.0;   % alpha_ms (D1)
MSparams(:,12) = 3.75;    % beta1_ms
MSparams(SIMPARAMS.net.MS.D2inds,12) = MSparams(SIMPARAMS.net.MS.D2inds,12) .* 0.0;   % beta1_ms (D2)
MSparams(:,13) = 0.156;   % beta2_ms
MSparams(SIMPARAMS.net.MS.D1inds,13) = MSparams(SIMPARAMS.net.MS.D1inds,13) .* 0.0;   % beta2_ms (D1)
MSparams(:,14) = 22.7;    % gDAms
SIMPAMS.MSparams(SIMPARAMS.net.MS.D2inds,14) = MSparams(SIMPARAMS.net.MS.D2inds,14) .* 0.0;   % gDA_ms (D2)

% -------------------------------------------------------------------------
% parameters for the FS neurons
% -------------------------------------------------------------------------
function FSparams = getFSparams(SIMPARAMS)
FSparams = ones(SIMPARAMS.net.FS.N,12); 

FSparams(:,1) = 80;       % C
FSparams(:,2) = -70;      % vr
FSparams(:,3) = -50;      % vt
FSparams(:,4) = 1;        % k
FSparams(:,5) = 0.2;      % a
FSparams(:,6) = 0.025;    % b
FSparams(:,7) = -60;      % c
FSparams(:,8) = 0.0;      % d
FSparams(:,9) = 25;       % vpeak
FSparams(:,10) = -55;     % vb
FSparams(:,11) = 0.1;     % eta
FSparams(:,12) = 0.625;   % epsilon

function network = LoadStriatumNetwork(fname)

load(fname);
if exist('SIMPARAMS');
    network = SIMPARAMS.net;
else
    network = net;
end


