%% script to assess basic firing stats of simulations...
clear all; close all

% load('../MSNonly.mat'); load('../MSNonly_SIMPARAMS.mat')
% load('../PW-like_MSNonly.mat'); load('../PW-like_MSNonly_SIMPARAMS.mat')
% load('../PW-like2_MSN_only.mat'); load('../PW-like2_MSN_only_SIMPARAMS.mat')
% load('../PW-like3_MSN_only.mat'); load('../PW-like3_MSN_only_SIMPARAMS.mat')
% load('../PW-likeSmall'); load('../PW-likeSmall_SIMPARAMS.mat')
% load('../PW-likeSmall2'); load('../PW-likeSmall2_SIMPARAMS.mat')
% load('../PW-likeSmall3.mat'); load('../PW-likeSmall3_SIMPARAMS.mat')
% load('../PW-likeSmall4.mat'); load('../PW-likeSmall4_SIMPARAMS.mat')
% load('../PW-likeSmall5.mat'); load('../PW-likeSmall5_SIMPARAMS.mat')
% load('../PW-likeSmall6.mat'); load('../PW-likeSmall6_SIMPARAMS.mat')
% load('../PW-likeSmall7.mat'); load('../PW-likeSmall7_SIMPARAMS.mat')
% load('../PW-likeSmall8.mat'); load('../PW-likeSmall8_SIMPARAMS.mat')
% load('../PW-likeSmall9.mat'); load('../PW-likeSmall9_SIMPARAMS.mat')

% load('../PW-likeVerySmall1.mat'); load('../PW-likeVerySmall1_SIMPARAMS.mat')
% load('../PW-likeVerySmall2.mat'); load('../PW-likeVerySmall2_SIMPARAMS.mat')
% load('../PW-likeVerySmall3.mat'); load('../PW-likeVerySmall3_SIMPARAMS.mat')
% load('../PW-like model results/PW-likeVerySmall4.mat'); load('../PW-like model results/PW-likeVerySmall4_SIMPARAMS.mat')

%% basic and swapped FSI dist models
% load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI1.mat'); load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI1_SIMPARAMS.mat')
% load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI3.mat'); load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI3_SIMPARAMS.mat')
% load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI5.mat'); load('../Basic and swapped FSI dist model results/250micron_model/BasicModel_FSI5_SIMPARAMS.mat')

% load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI1.mat'); load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI1_SIMPARAMS.mat');
% load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI3.mat'); load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI3_SIMPARAMS.mat');
% load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI5.mat'); load('../Basic and swapped FSI dist model results/250micron_model/SwappedFSIdists_FSI5_SIMPARAMS.mat');

% 500 micron cubes
% load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI1.mat'); load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI1_SIMPARAMS.mat')
% load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI3.mat'); load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI3_SIMPARAMS.mat')
% load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI5.mat'); load('../Basic and swapped FSI dist model results/500micron_model/BasicModel_FSI5_SIMPARAMS.mat')

% load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI1.mat'); load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI1_SIMPARAMS.mat');
% load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI3.mat'); load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI3_SIMPARAMS.mat');
load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI5.mat'); load('../Basic and swapped FSI dist model results/500micron_model/SwappedFSIdists_FSI5_SIMPARAMS.mat');

%load('../Sphere_100micron_radius.mat');
% load('../Sphere_150micron_radius.mat');
% load('../Sphere_100micron_radius_no_bghz.mat');
% load('../Sphere_100micron_radius_1pt5_bghz.mat');
% load('../Sphere_150micron_radius_1pt5_bghz.mat');
% load('../Sphere_200micron_radius_1pt5_bghz.mat');
% load('../Sphere_150micron_radius_MSNFSIs_0_bghz.mat');
% load('../Sphere_250micron_radius_MSNFSIs_0_bghz.mat');
%% effects on BG input responding MSNs...
% load('../Sphere_no_input_MSNFSIs_1pt9_bghz.mat');
% load('../Sphere_150micron_radius_MSNFSIs_1pt9_bghz.mat');
% load('../Sphere_250micron_radius_MSNFSIs_1pt9_bghz.mat');

MSspks = out.STms; 
MSspks(:,1) = MSspks(:,1)+1; % change from zero-base to 1-base index 
MSidxs = unique(MSspks(:,1));
Nmsidxs = numel(MSidxs);
Nms = SIMPARAMS.net.MS.N;

FSspks = out.STfs;
FSspks(:,1) = FSspks(:,1)+1; % change from zero-base to 1-base index 
FSidxs = unique(FSspks(:,1));
Nfsidxs = numel(FSidxs);
Nfs = SIMPARAMS.net.FS.N;

% simulation parameters
simT = SIMPARAMS.sim.tfinal; % in ms
T = simT * 1e-3;    % in seconds

% get stats for all MSNs
MSrate = zeros(Nms,1); MS_CVisi = zeros(Nms,1);
for j = 1:Nmsidxs
    currix = find(MSspks(:,1) == MSidxs(j));
    % basic rate
    MSrate(MSidxs(j)) = numel(currix) / T;
    
    % ISIs
    ts = MSspks(currix,2)*1e-3; % in seconds
    ISIs = diff(ts);
    MS_CVisi(MSidxs(j)) = std(ISIs)/mean(ISIs);
end
MS_CVisi(isnan(MS_CVisi)) = 0;

% get stats for all FSIs
FSrate = zeros(Nfs,1); FS_CVisi = zeros(Nfs,1);
for j = 1:Nfsidxs
    currix = find(FSspks(:,1) == FSidxs(j));
    % basic rate
    FSrate(FSidxs(j)) = numel(currix) / T;
    
    % ISIs
    ts = FSspks(currix,2)*1e-3; % in seconds
    ISIs = diff(ts);
    FS_CVisi(FSidxs(j)) = std(ISIs)/mean(ISIs);
end
FS_CVisi(isnan(FS_CVisi)) = 0;

%%%% plot empirical CDFs... 
%%% NOTE that stats structure returned by cdfplot has all info you need in
%%% it for all neurons
% with addition of:
MSrate_IQR = iqr(MSrate);
FSrate_IQR = iqr(FSrate);
MSCV_IQR = iqr(MS_CVisi);
FSCV_IQR = iqr(FS_CVisi);

figure(10); clf
[hM,MSstats] = cdfplot(MSrate); hold on
[hF,FSstats] = cdfplot(FSrate); set(hF,'Color',[1 0 0])
phF = get(hF,'Parent'); set(phF,'XScale','log')
title('ECDF plots for firing rates')
xlabel('x=firing rate (spikes/s)')

[MSrt_ecdf,MSrt_x] = ecdf(MSrate);
[FSrt_ecdf,FSrt_x] = ecdf(FSrate);

figure(101); clf
subplot(211), hist(MSrate,logspace(-1,2,25))
%bar(logspace(-2,2,25),N)
subplot(212), hist(FSrate,logspace(-1,2,25))

figure(11); clf
[hMCV,MSCVstats] = cdfplot(MS_CVisi); hold on
[hFCV,FSCVstats] = cdfplot(FS_CVisi(~isnan(FS_CVisi))); set(hFCV,'Color',[1 0 0])
title('ECDF plots for ISI CV')
xlabel('x=ISI CV')

[MSCV_ecdf,MSCV_x] = ecdf(MS_CVisi);
[FSCV_ecdf,FSCV_x] = ecdf(FS_CVisi);


% array of data for Sigmaplot sheet
data_vec = [MSstats.median MSrate_IQR FSstats.median FSrate_IQR MSCVstats.median MSCV_IQR FSCVstats.median FSCV_IQR]'
% data_vec_means = [MSstats.mean MSstats.std FSstats.mean FSstats.std MSCVstats.median MSCV_IQR FSCVstats.median FSCV_IQR]'

% characteristics of active cells!
FS_med_nz_rt = median(FSrate(FSrate>0))
MS_med_nz_rt = median(MSrate(MSrate>0))
FS_med_nz_CV = median(FS_CVisi(FS_CVisi>0))
MS_med_nz_CV = median(MS_CVisi(MS_CVisi>0))

% plot raster in CV or rate order....
%[x I] = sort(MS_CVisi); 
[x I] = sort(MSrate);
tempix = zeros(numel(MSspks(:,1)),1);
for i =1:Nms
    ixs = find(MSspks(:,1) == I(i));
    tempix(ixs) = i;    
end
h = raster_plot(tempix,MSspks(:,2));

[x I] = sort(FS_CVisi); 
% [x I] = sort(FSrate);
FStempix = zeros(numel(FSspks(:,1)),1);
for i =1:Nfs
    ixs = find(FSspks(:,1) == I(i));
    FStempix(ixs) = i;    
end
FStempix = FSspks(:,1);

t = [2000 3000]; % section to look at in ms
ixs = find(FSspks(:,2) >= t(1) & FSspks(:,2)<= t(2));
h2 = raster_plot(FStempix(ixs,1),FSspks(ixs,2)); c2 = get(h2,'Children');
axis([t(1) t(2) 0 Nfs+1]); set(c2,'XTick',[],'YTick',[]); xlabel(''); ylabel('');

% pD1 = mean(nD1 ./ (nD1+nD2));

% but we want to analyse just *active cells* ?
% MSrate_IQRb = iqr(MSrate(MSrate>0));
% FSrate_IQRb = iqr(FSrate(FSrate>0));
% MSCV_IQRb = iqr(MS_CVisi(MS_CVisi > 0));
% FSCV_IQRb = iqr(FS_CVisi(FS_CVisi > 0));
% MSmedian = median(MSrate(MSrate>0));
% FSmedian = median(FSrate(FSrate>0));
% MSCV_median = median(MS_CVisi(MS_CVisi > 0));
% FSCV_median = median(FS_CVisi(FS_CVisi > 0));
% 
% data_vecb = [MSmedian MSrate_IQRb FSmedian FSrate_IQRb MSCV_median MSCV_IQRb FSCV_median FSCV_IQRb]';
save MSrt_x.txt MSrt_x -ascii
save MSrt_ecdf.txt MSrt_ecdf -ascii
save FSrt_x.txt FSrt_x -ascii
save FSrt_ecdf.txt FSrt_ecdf -ascii
save MSCV_x.txt MSCV_x -ascii
save MSCV_ecdf.txt MSCV_ecdf -ascii
save FSCV_x.txt FSCV_x -ascii
save FSCV_ecdf.txt FSCV_ecdf -ascii




