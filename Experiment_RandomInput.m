function out = Experiment_RandomInput(DA, fname)
if nargin == 0
    DA = 0.0;
    fname = {'RandomInput'};
end

% set the model parameters
SIMPARAMS = StriatumNetworkParameters;

% name for the log file
SIMPARAMS.sim.logfname = [char(fname) '.log'];

% -------------------------------------------------------------------------
% set the DA level
SIMPARAMS.physiology.DA = DA; 

% these will need re-scaling to account for removal of 1/tau_s if using RK2_B engine 
SIMPARAMS.net.Cctms_w = ones(length(SIMPARAMS.net.Cctms),1) .* 6.1;   
SIMPARAMS.net.Cctfs_w = ones(length(SIMPARAMS.net.Cctfs),1) .* 6.1; 
SIMPARAMS.net.Cmsms_w = ones(length(SIMPARAMS.net.Cmsms),1) .* 4.36; 
SIMPARAMS.net.Cfsms_w = ones(length(SIMPARAMS.net.Cfsms),1) .* (4.36 * 5); 
SIMPARAMS.net.Cfsfs_w = ones(length(SIMPARAMS.net.Cfsfs),1) .* (4.36 * 5); % same as FS-MSN input....
SIMPARAMS.net.Cgapfs_w = ones(length(SIMPARAMS.net.Cgapfs_w),1).* (150/5); % too strong?

% set the input parameters - overrides settings in StriatumNetworkParameters 
SIMPARAMS.sim.tfinal = 10000; % msec
% SIMPARAMS.input.CTX.r_MSSEG = ones(SIMPARAMS.net.MS.N,1) .* 1.9; % Hz
% SIMPARAMS.input.CTX.N_MSSEG = int32(ones(SIMPARAMS.net.MS.N,1) .* 250);
% SIMPARAMS.input.CTX.r_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* 1.9; % Hz
% SIMPARAMS.input.CTX.N_FSSEG = int32(ones(SIMPARAMS.net.FS.N,1) .* 250);
% SIMPARAMS.input.CTX.alpha_FSSEG = ones(SIMPARAMS.net.FS.N,1) .* 0.0;

save([char(fname) '_SIMPARAMS'],'SIMPARAMS');
% -------------------------------------------------------------------------
% Run the simulation
tic
out = RunSimulation(SIMPARAMS);
toc

% -------------------------------------------------------------------------
% Save the results to disc
save(char(fname), 'out');

% -------------------------------------------------------------------------
% plot the results
close all
figure(1); clf; plot(out.STms(:,2), out.STms(:,1), '.')
figure(2); clf; plot(out.STfs(:,2), out.STfs(:,1), '.')
figure(3); clf; plot(out.RecordChan_MS(:,1:25))
figure(4); clf; plot(out.RecordChan_MS(:,26:end))
