% Script containing analyses for Figure 5b of "Dynamic hidden
% states underlying working-memory-guided behavior", Nature Neuroscience,
% 2017

% Input is preprocessed (see article) and baselined (-200 to 0 ms relative
% to stimulus onset) EEG data.

close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states'; % path to main folder containing all data, functions, toolboxes, scripts
addpath(genpath(main_dir))
dat_dir=[main_dir '/Data'];

angspace=(-pi:pi/6:pi)'; % angular space for tuning curve (in radians)
angspace(end)=[];
bin_width=pi/6; % width of each angle bin of the tuning curve (in radians)

s_factor=8; % determines the SD of the smoothing kernel for the decoding time-courses,
            % depends on resolution of data (500 hz = 2ms, 2ms times 8 is thus an SD of 16 ms)
%%
for sub=1:19
    
    fprintf(['Doing ' num2str(sub) '\n'])
    
    load(fullfile(dat_dir,['Dynamic_hidden_states_exp2_' num2str(sub) '.mat']));
    
    EEG_dat1=exp2_data.EEG_mem_items_sess1;
    Results1=exp2_data.Results_sess1;
    
    EEG_dat2=exp2_data.EEG_mem_items_sess2;
    Results2=exp2_data.Results_sess2;
    
    time=EEG_dat1.time;
    
    clear exp2_data
    
    incl1=not(ismember(1:size(EEG_dat1.trial,1),EEG_dat1.bad_trials))'; % logical array of trials to be included
    data1 = EEG_dat1.trial(incl1,:,:);
    
    incl2=not(ismember(1:size(EEG_dat2.trial,1),EEG_dat2.bad_trials))'; % logical array of trials to be included
    data2 = EEG_dat2.trial(incl2,:,:); 
    
    clear EEG_dat1 EEG_dat2
    
    data1= bsxfun(@minus, data1, mean(data1,2)); % mean center voltage across channels
    mem_angles1=Results1(incl1,1:2)*2; % extract memory item angles and rescale 
    
    data2= bsxfun(@minus, data2, mean(data2,2)); % mean center voltage across channels
    mem_angles2=Results2(incl2,1:2)*2; % extract memory item angles and rescale

    dec_early1 = mahalTune_func_cross_temp(data1,mem_angles1(:,1),angspace,bin_width); % cross-temp decoding of early-tested item in session 1
    dec_late1 = mahalTune_func_cross_temp(data1,mem_angles1(:,2),angspace,bin_width); % cross-temp decoding of late-tested item in session 1
    
    dec_early2 = mahalTune_func_cross_temp(data2,mem_angles2(:,1),angspace,bin_width); % cross-temp decoding of early-tested item in session 2
    dec_late2 = mahalTune_func_cross_temp(data2,mem_angles2(:,2),angspace,bin_width); % cross-temp decoding of late-tested item in session 2
    
    dec_mem_early(sub,:,:)   = imgaussfilt(((dec_early1+dec_early2)/2 + permute((dec_early1+dec_early2)/2,[2 1]))./2,s_factor); % average over both sessions and make symmetrical
    dec_mem_late(sub,:,:)   = imgaussfilt(((dec_late1+dec_late2)/2 + permute((dec_late1+dec_late2)/2,[2 1]))./2,s_factor);
end

%% prepare data and test for dynamic coding
% test if on-diagonal decoding is significantly higher than
% off-diagonal decoding
%% early item
data=dec_mem_early(:,time>0,time>0);

nsubs  = size(data,1);
ntimes = size(data,2);
xdiag  = nan(nsubs,ntimes,ntimes);
ydiag  = nan(nsubs,ntimes,ntimes);

for itime = 1:ntimes
    %get on-diagonal decoding value
    currdat =  data(:,itime,itime);
    %create a matrix with the on-diagonal element along all points in a row
    %(i.e. time axis 1)
    xdiag(:,itime,:) = repmat(currdat,[1 ntimes]);
    %create a matrix with the on-diagonal element along all points in a
    %column (i.e. time axis 2)
    ydiag(:,:,itime) = repmat(currdat,[1 ntimes]);
end

%subtract on-diagonal from off-diagonal values along both time axes
xdiff = data - xdiag;
ydiff = data - ydiag;

%because matrix is symmetrical, only keep the upper triangular part of it
for sub = 1:nsubs
    xdiff(sub,:,:) = triu(squeeze(xdiff(sub,:,:)));
    ydiff(sub,:,:) = triu(squeeze(ydiff(sub,:,:)));
end
 
% do cluster-corrected test late-tested item
pclust_early = ClusterCorrectionConjunction(xdiff,ydiff,5000,0.05,0.05,false); 
h_clust_early=squeeze(pclust_early<0.05);

%% late item
data=dec_mem_late(:,time>0,time>0);

nsubs  = size(data,1);
ntimes = size(data,2);
xdiag  = nan(nsubs,ntimes,ntimes);
ydiag  = nan(nsubs,ntimes,ntimes);

for itime = 1:ntimes
    %get on-diagonal decoding value
    currdat =  data(:,itime,itime);
    %create a matrix with the on-diagonal element along all points in a row
    %(i.e. time axis 1)
    xdiag(:,itime,:) = repmat(currdat,[1 ntimes]);
    %create a matrix with the on-diagonal element along all points in a
    %column (i.e. time axis 2)
    ydiag(:,:,itime) = repmat(currdat,[1 ntimes]);
end

%subtract on-diagonal from off-diagonal values along both time axes
xdiff = data - xdiag;
ydiff = data - ydiag;

%because matrix is symmetrical, only keep the upper triangular part of it
for sub = 1:nsubs
    xdiff(sub,:,:) = triu(squeeze(xdiff(sub,:,:)));
    ydiff(sub,:,:) = triu(squeeze(ydiff(sub,:,:)));
end
 
% do cluster-corrected test late-tested item
pclust_late = ClusterCorrectionConjunction(xdiff,ydiff,5000,0.05,0.05,false); 
h_clust_late=squeeze(pclust_late<0.05);

%% plot dynamics and significant dynamic clusters

% early-tested item
B = bwboundaries(h_clust_early);
fhandle=figure(1);
title('Figure 5B, left')
imagesc(time,time,squeeze(mean(dec_mem_early,1)),[-0.001 0.001]); axis xy
hold all
for k=1:length(B)
    boundary=B{k};
    plot(((boundary(:,2)+sum(time<=0))/500-.100),((boundary(:,1)+sum(time<=0))/500-.100),'k')
end
colormap('jet')
xlim([-0.1 1.2]);ylim([-0.1 1.2])
pbaspect([1,1,1]);set(gca,'TickDir','out');
xlabel('Time (s) - relative to memory items')
ylabel('Time (s) - relative to memory items')

% late-tested item
B = bwboundaries(h_clust_late);
fhandle=figure(2);
title('Figure 5B, middle')
imagesc(time,time,squeeze(mean(dec_mem_late,1)),[-0.001 0.001]); axis xy
hold all
for k=1:length(B)
    boundary=B{k};
    plot(((boundary(:,2)+sum(time<=0))/500-.100),((boundary(:,1)+sum(time<=0))/500-.100),'k')
end
colormap('jet')
xlim([-0.1 1.2]);ylim([-0.1 1.2])
pbaspect([1,1,1]);set(gca,'TickDir','out');
xlabel('Time (s) - relative to memory items')
ylabel('Time (s) - relative to memory items')

%% cluster-corrected stats on the difference between items
dec_diff=dec_mem_early-dec_mem_late; % difference between early and late tested item decoding
[datobs, datrnd] = cluster_test_helper(permute(dec_diff(:,time>0,time>0),[2 3 1]), 5000);
[h, p_mem_diff, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

% plot the difference and corresponding significant clusters
B = bwboundaries(h);
fhandle=figure(3);
title('Figure 5B, right')
imagesc(time,time,squeeze(mean(dec_diff,1)),[-0.001 0.001]); axis xy
hold all
for k=1:length(B)
    boundary=B{k};
    plot(((boundary(:,2)+sum(time<=0))/500-.100),((boundary(:,1)+sum(time<=0))/500-.100),'k')
end
colormap('jet')
xlim([-0.1 1.2]);ylim([-0.1 1.2])
pbaspect([1,1,1]);set(gca,'TickDir','out');
xlabel('Time (s) - relative to memory items')
ylabel('Time (s) - relative to memory items')

