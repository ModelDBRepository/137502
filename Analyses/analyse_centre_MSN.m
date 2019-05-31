%% script to assess basic firing stats of simulations...
clear all; close all;

width = 50;
% radius = [0 100 200 300 400];
radius = 0:50:450;

iqrMSNrate = zeros(numel(radius),2); iqrMSNCVisi = zeros(numel(radius),2);
iqrFSrate = zeros(numel(radius),2); iqrFSCVisi = zeros(numel(radius),2);
for loop = 1:numel(radius)
    % fname = ['../Shell input results/500micron_cube/Shell_width' num2str(width) '_innerradius_' num2str(radius(loop)) '.mat'];
    % fname = ['../Shell input results/1mm cube/3%_FSIs/Shell_width' num2str(width) '_innerradius_' num2str(radius(loop)) '.mat']; % 1mm cubed....
    fname = ['../Shell input results/1mm cube/1%_FSIs/Shell_width' num2str(width) '_innerradius_' num2str(radius(loop)) '.mat']; % 1mm cubed....

    load(fname);

    MSspks = out.STms; 
    MSspks(:,1) = MSspks(:,1)+1; % change from zero-base to 1-base index 
    Nmsidxs(loop) = numel(SIMPARAMS.input.shell.MSids);   
    Nms_all = SIMPARAMS.net.MS.N;

    FSspks = out.STfs;
    FSspks(:,1) = FSspks(:,1)+1; % change from zero-base to 1-base index 
    Nfsidxs(loop) = numel(SIMPARAMS.input.shell.FSids);
    Nfs_all = SIMPARAMS.net.FS.N;

    % simulation parameters
    simT = SIMPARAMS.sim.tfinal; % in ms
    T = simT * 1e-3;    % in seconds

    % get stats for all MSNs other than centre
    MSrate = zeros(Nmsidxs(loop),1); MS_CVisi = zeros(Nmsidxs(loop),1);
    for j = 1:Nmsidxs(loop)
        if j ~= SIMPARAMS.input.shell.MScentre
            currix = find(MSspks(:,1) == SIMPARAMS.input.shell.MSids(j));
            % basic rate
            MSrate(j) = numel(currix) / T;

            % ISIs
            ts = MSspks(currix,2)*1e-3; % in seconds
            ISIs = diff(ts);
            MS_CVisi(j) = std(ISIs)/mean(ISIs);
        end
    end
    MS_CVisi(isnan(MS_CVisi)) = 0;
    medianMSNrate(loop) = median(MSrate); medianMSNCVisi(loop) = median(MS_CVisi); meanMSNrate(loop) = mean(MSrate); 
    iqrMSNrate(loop,:) = prctile(MSrate,[25,75]); iqrMSNCVisi(loop,:) = prctile(MS_CVisi,[25,75]);
    
    
    % get stats for all FSIs
    FSrate = zeros(Nfsidxs(loop),1); FS_CVisi = zeros(Nfsidxs(loop),1);
    for j = 1:Nfsidxs(loop)
        currix = find(FSspks(:,1) == SIMPARAMS.input.shell.FSids(j));
        % basic rate
        FSrate(j) = numel(currix) / T;

        % ISIs
        ts = FSspks(currix,2)*1e-3; % in seconds
        ISIs = diff(ts);
        FS_CVisi(j) = std(ISIs)/mean(ISIs);
    end
    FS_CVisi(isnan(FS_CVisi)) = 0;
    medianFSrate(loop) = median(FSrate); medianFSCVisi(loop) = median(FS_CVisi);  meanFSrate(loop) = mean(FSrate); 
    iqrFSrate(loop,:) = prctile(FSrate,[25,75]); iqrFSCVisi(loop,:) = prctile(FS_CVisi,[25,75]);

    %%% centre MSN
    currix = find(MSspks(:,1) == SIMPARAMS.input.shell.MScentre);
    centre_rate(loop) = numel(currix) / T;
    ts = MSspks(currix,2)*1e-3; % in seconds
    ISIs = diff(ts);
    centre_CVisi(loop) = std(ISIs)/mean(ISIs);

    %%% number of connected neurons to centre within shell....
    MSNcnctd = [];
    for j = 1:Nmsidxs(loop)
        % find all MSNs that MSN contacts...
        thisID = SIMPARAMS.input.shell.MSids(j);
        tgts = SIMPARAMS.net.Cmsms(SIMPARAMS.net.Cmsms_b(thisID)+1:SIMPARAMS.net.Cmsms_b(thisID+1))+1;  % add 1 to index cos is 0-base
        blnCnct = find(tgts == SIMPARAMS.input.shell.MScentre);
        if blnCnct MSNcnctd = [MSNcnctd; thisID]; end
    end
    
    FSIcnctd = [];
    for j = 1:Nfsidxs(loop)
        % find all MSNs that MSN contacts...
        thisID = SIMPARAMS.input.shell.FSids(j);
        tgts = SIMPARAMS.net.Cfsms(SIMPARAMS.net.Cfsms_b(thisID)+1:SIMPARAMS.net.Cfsms_b(thisID+1))+1;  % add 1 to index cos is 0-base
        blnCnct = find(tgts == SIMPARAMS.input.shell.MScentre);
        if blnCnct FSIcnctd = [FSIcnctd; thisID]; end
    end

    nMSNcnctd(loop) = numel(MSNcnctd); nFSIcnctd(loop) = numel(FSIcnctd);
    pMSNcnctd(loop) = nMSNcnctd(loop) / Nmsidxs(loop); pFSIcnctd(loop) = nFSIcnctd(loop) / Nfsidxs(loop);
end

figure(1); clf
subplot(311); plot(radius,Nmsidxs,'+-'); hold on; plot(radius,Nfsidxs,'r+-');
xlabel('Inner sphere radius (microns)'); ylabel('Number of neurons in shell');
subplot(312); plot(radius,nMSNcnctd,'+-'); hold on; plot(radius,nFSIcnctd,'r+-');
xlabel('Inner sphere radius (microns)'); ylabel('Number of neurons in shell projecting to centre MSN');
subplot(313); plot(radius,pMSNcnctd,'+-'); hold on; plot(radius,pFSIcnctd,'r+-');
xlabel('Inner sphere radius (microns)'); ylabel('Proportion of neurons in shell projecting to centre MSN');

figure(2); clf
subplot(211); plot(radius,medianMSNrate,'+-'); hold on; plot(radius,iqrMSNrate(:,1),'+:');  plot(radius,iqrMSNrate(:,2),'+:');
plot(radius,medianFSrate,'r+-'); plot(radius,iqrFSrate(:,1),'r+:');  plot(radius,iqrFSrate(:,2),'r+:'); 
xlabel('Inner sphere radius (microns)'); ylabel('Firing rate of neurons projecting to centre MSN (spikes/s)');
subplot(212); plot(radius,medianMSNCVisi,'+-'); hold on; plot(radius,iqrMSNCVisi(:,1),'+:');  plot(radius,iqrMSNCVisi(:,2),'+:');
plot(radius,medianFSCVisi,'r+-'); plot(radius,iqrFSCVisi(:,1),'r+:');  plot(radius,iqrFSCVisi(:,2),'r+:'); 
xlabel('Inner sphere radius (microns)'); ylabel('ISI CV of neurons projecting to centre MSN');

figure(3); clf
subplot(211); plot(radius,centre_rate,'+-');
xlabel('Inner sphere radius (microns)'); ylabel('centre MSN rate (spikes/s)');
subplot(212); plot(radius,centre_CVisi,'+-'); 
xlabel('Inner sphere radius (microns)'); ylabel('centre MSN ISI CV');


save 1prct_FSI_centre_MSN

