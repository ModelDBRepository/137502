function out = Experiment_ImpactOnCentralMSN(DA,width,radius)
if nargin == 0
    width = 50;
    radius = 100;
    DA = 0.0;
    fname = {'Expt_InputsToCentreMSN'};
else
    fname = {['Shell_width' num2str(width) '_innerradius_' num2str(radius)]}; 
end
% set the model parameters
% network_is_at = 'C:\Users\Mark\Simulations\Striatum models\Simulations for PLoS striatal network paper\Sphere_150micron_radius_MSNFSIs_0_bghz.mat';
% network_is_at = 'C:\Simulations\Striatum models\Simulations for PLoS striatal network paper\CreateNetwork\Striatum_network_1000-1000-1000_num_1_at_734300.5522.mat'; % 3% FSI network
network_is_at = 'C:\Simulations\Striatum models\Simulations for PLoS striatal network paper\CreateNetwork\Striatum_network_1000-1000-1000_num_1_at_734299.3923.mat'; % 1% FSI network
SIMPARAMS = StriatumNetworkParameters(network_is_at);   % is vital we use the same network each time!!!
% SIMPARAMS = StriatumNetworkParameters;

% name for the log file
SIMPARAMS.sim.logfname = [char(fname) '.log'];

% -------------------------------------------------------------------------
% set the DA level
SIMPARAMS.physiology.DA = DA; 

% set all the GABA the weights to 0
SIMPARAMS.net.Cctms_w = ones(length(SIMPARAMS.net.Cctms),1) .* 6.1;
SIMPARAMS.net.Cctfs_w = ones(length(SIMPARAMS.net.Cctfs),1) .* 6.1;
SIMPARAMS.net.Cmsms_w = ones(length(SIMPARAMS.net.Cmsms),1) .* 4.36;
SIMPARAMS.net.Cfsms_w = ones(length(SIMPARAMS.net.Cfsms),1) .* (4.36 * 5);
SIMPARAMS.net.Cfsfs_w = ones(length(SIMPARAMS.net.Cfsfs),1) .* (4.36 * 5);
SIMPARAMS.net.Cgapfs_w = ones(length(SIMPARAMS.net.Cgapfs_w),1).* (150/5); 

% set the backgroung input level
bgHz = 0;   % normally 1.9; only interested in effect of shell input here... 
SIMPARAMS.input.CTX.r_MSSEG = ones(SIMPARAMS.net.MS.N,1) .* bgHz;
SIMPARAMS.input.CTX.N_MSSEG = int32(ones(SIMPARAMS.net.MS.N,1) .* 250);
SIMPARAMS.input.CTX.alpha_MSSEG = ones(SIMPARAMS.net.MS.N,1) .* 0.0;
SIMPARAMS.input.CTX.r_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* bgHz;
SIMPARAMS.input.CTX.N_FSSEG = int32(ones(SIMPARAMS.net.FS.N,1) .* 250);
SIMPARAMS.input.CTX.alpha_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* 0.0;

% standard selection sims....
SIMPARAMS.sim.tfinal = 4000;
hz_input = 5; 
% SIMPARAMS.input.Selection.Pt =  int32([500 550 750 800] ./ SIMPARAMS.sim.dt)'; % time of change in msec, converted into iteration
% SIMPARAMS.input.Selection.Pch = int32([1 2 1 2])'; % channel
% SIMPARAMS.input.Selection.Phz = [4 6 bgHz bgHz]'; % rate of sp

SIMPARAMS.input.shell.inner_sphere_radius = radius; %n_inputs = 200; nFS_inputs = 20; % assuming minimum sphere diameter is 300 microns, using 3% FSIs....
SIMPARAMS.input.shell.shell_width = width;
mid = SIMPARAMS.net.PhysicalDimensions/2; % centre of cubic network  
MSdist_from_mid = sqrt(sum((SIMPARAMS.net.MS.Position - repmat(mid,SIMPARAMS.net.MS.N,1))'.^2));
FSdist_from_mid = sqrt(sum((SIMPARAMS.net.FS.Position - repmat(mid,SIMPARAMS.net.FS.N,1))'.^2));

SIMPARAMS.input.shell.MScentre = find(MSdist_from_mid == min(MSdist_from_mid));
SIMPARAMS.input.shell.MSids = find(MSdist_from_mid >= SIMPARAMS.input.shell.inner_sphere_radius & MSdist_from_mid <= SIMPARAMS.input.shell.inner_sphere_radius + SIMPARAMS.input.shell.shell_width);   % all MSNs within this shell
SIMPARAMS.input.shell.FSids = find(FSdist_from_mid >= SIMPARAMS.input.shell.inner_sphere_radius & FSdist_from_mid <= SIMPARAMS.input.shell.inner_sphere_radius + SIMPARAMS.input.shell.shell_width);   % all MSNs within this shell

SIMPARAMS.input.CTX.r_MSSEG(SIMPARAMS.input.shell.MScentre) = hz_input; % centre MSN....
SIMPARAMS.input.CTX.r_MSSEG(SIMPARAMS.input.shell.MSids) = hz_input; % all MSNs get input in this shell
SIMPARAMS.input.CTX.r_FSSEG(SIMPARAMS.input.shell.FSids) = hz_input; % all FSIs get input in this shell

% keyboard

% in case we want to randomly assign inputs....
%tmp = randperm(numel(MSids));
%tmp = randperm(numel(FSids));
% SIMPARAMS.input.CTX.r_MSSEG(MSids(tmp(1:n_inputs))) = hz_input; % randomly assign which MSNs get input in this sphere...
% SIMPARAMS.input.CTX.r_FSSEG(FSids(tmp(1:nFS_inputs))) = hz_input; % randomly assign which MSNs get input in this sphere...

% frac_MSNs = n_inputs ./ ((4/3 * pi * (sphere_width/1000)^3) * SIMPARAMS.net.CellsPerMillimeterCube); 
% nFS_inputs = round(SIMPARAMS.net.FS.N * frac_MSNs);


% Phasic input to MSNs only
SIMPARAMS.net.CHAN1_MS = int32(0:(round(SIMPARAMS.net.MS.N / 2)-1))';
SIMPARAMS.net.CHAN1_FS = int32(0:(round(SIMPARAMS.net.FS.N / 2)-1))' .* 0.0;  
SIMPARAMS.net.CHAN2_MS = int32(round(SIMPARAMS.net.MS.N / 2):(SIMPARAMS.net.MS.N-1))';
SIMPARAMS.net.CHAN2_FS = int32(round(SIMPARAMS.net.FS.N / 2):(SIMPARAMS.net.FS.N-1))' .* 0.0;



% -------------------------------------------------------------------------
% Run the simulation
tic
out = RunSimulation(SIMPARAMS);
toc

% -------------------------------------------------------------------------
% Save the results to disc
save(char(fname), 'out', 'SIMPARAMS') ;

% -------------------------------------------------------------------------
% Sort and plot the results
close all

% -------------------------------------------------------------------------
close all
figure(1); clf; plot(out.STms(:,2), out.STms(:,1), '.')
figure(2); clf; plot(out.STfs(:,2), out.STfs(:,1), '.')

% sort out the spikes according to channel and plot the results

