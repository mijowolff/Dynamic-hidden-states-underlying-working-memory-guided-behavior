% Script containing analyses for Figures 2f, 3a, and 3b of "Dynamic hidden
% states underlying working-memory-guided behavior", Nature Neuroscience,
% 2017

% Input is preprocessed (see article) and baselined (-200 to 0 ms relative
% to stimulus onset) EEG data.

close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states'; % path to main folder containing all data, functions/toolboxes, scripts...
addpath(genpath(main_dir))
dat_dir=[main_dir '/Data'];

angspace=(-pi:pi/6:pi)'; % angular space for tuning curve (in radians)
angspace(end)=[];
bin_width=pi/6; % width of each angle bin of the tuning curve (in radians)

s_factor=8; % determines the SD of the smoothing kernel for the decoding time-courses,
% depends on resolution of data (500 hz = 2ms, 2ms times 8 is thus an SD of 16 ms)

%%
for sub=1:30
    fprintf(['Doing ' num2str(sub) '\n'])
    
    fname_in=['Dynamic_hidden_states_exp1_' num2str(sub) '.mat'];
    
    load(fullfile(dat_dir,fname_in));
    
    EEG_dat=exp1_data.EEG_impulse;
    Results=exp1_data.Results;
    time=EEG_dat.time;
    
    clear exp1_data
    
    incl=not(ismember(1:size(EEG_dat.trial,1),EEG_dat.bad_trials))'; % logical array of trials to be included
    
    data=EEG_dat.trial(incl,:,:);
    
    data= bsxfun(@minus, data, mean(data,2)); % mean center voltage across channels to normalize
    
    mem_angles=Results(incl,1:2)*2; % extract memory item angles and rescale
    cue=Results(incl,3);
    acc=Results(incl,5);
    probe_rotation=Results(incl,4);
    
    cued_mem_left=mem_angles(cue==1,1); %orientations of cued items on the left
    cued_mem_right=mem_angles(cue==2,2); %orientations of cued items on the right
    
    uncued_mem_left=mem_angles(cue==2,1); %orientations of uncued items on the left
    uncued_mem_right=mem_angles(cue==1,2); %orientations of uncued items on the right
    
    dec_cued_left = mahalTune_func(data(cue==1,:,:),cued_mem_left,angspace,bin_width); % decode cued items on the left
    dec_cued_right = mahalTune_func(data(cue==2,:,:),cued_mem_right,angspace,bin_width); % decode cued items on the right
    
    dec_imp_cued(sub,:)=imgaussfilt(mean(cat(1,dec_cued_left,dec_cued_right),1),s_factor); % average over trials and cued locations and smooth over time
    
    dec_uncued_left = mahalTune_func(double(data(cue==2,:,:)),uncued_mem_left,angspace,bin_width); % decode uncued items on the left
    dec_uncued_right = mahalTune_func(double(data(cue==1,:,:)),uncued_mem_right,angspace,bin_width); % decode uncued items on the right
    
    dec_imp_uncued(sub,:)=imgaussfilt(mean(cat(1,dec_uncued_left,dec_uncued_right),1),s_factor);
    
    clear EEG_dat
    %% the following section extracts high and low decoding trials for figure 3
    
    acc=cat(1,acc(cue==1,1),acc(cue==2)); % put the trials in the same order as the decoding output
    probe_rotation=cat(1,probe_rotation(cue==1,1),probe_rotation(cue==2,1));
    
    dec_cued_both_av=mean(cat(1,dec_cued_left(:,time>.1),dec_cued_right(:,time>.1)),2); %average decoding time-course (from 100 ms to end of trial)
    dec_uncued_both_av=mean(cat(1,dec_uncued_right(:,time>.1),dec_uncued_left(:,time>.1)),2);
    
    high_dec_cued_ind=dec_cued_both_av>median(dec_cued_both_av); %median split high...
    low_dec_cued_ind=dec_cued_both_av<median(dec_cued_both_av); %... and low decoding trials
    
    high_dec_uncued_ind=dec_uncued_both_av>median(dec_uncued_both_av);
    low_dec_uncued_ind=dec_uncued_both_av<median(dec_uncued_both_av);
    
    acc_diff_high_low_cued(sub,1)=mean(acc(high_dec_cued_ind,1),1)-mean(acc(low_dec_cued_ind,1),1); % accuracy difference between high and low decoding trials
    acc_diff_high_low_uncued(sub,1)=mean(acc(high_dec_uncued_ind,1),1)-mean(acc(low_dec_uncued_ind,1),1);
    
    resp = acc; resp(resp==0)=-1;resp = resp.* sign(probe_rotation);resp(resp==-1)=0; %convert accuracy to cw (1) or acw (0) response for behavioural modelling later on
    
    resp_high_cued=resp(high_dec_cued_ind,1);
    resp_high_uncued=resp(high_dec_uncued_ind,1);
    resp_low_cued=resp(low_dec_cued_ind,1);
    resp_low_uncued=resp(low_dec_uncued_ind,1);
    
    % prepare for modelling
    c=1;
    for p=unique(probe_rotation)'
        num_imp_cued_high(sub,c)=sum(resp_high_cued(probe_rotation(high_dec_cued_ind,1)==p,1),1); % number of clockwise responses for each probe-rotation condition separately
        num_imp_cued_low(sub,c)=sum(resp_low_cued(probe_rotation(low_dec_cued_ind,1)==p,1),1);
        tot_imp_cued_high(sub,c)=sum(probe_rotation(high_dec_cued_ind,1)==p); % total number of specific probe-rotaion condition
        tot_imp_cued_low(sub,c)=sum(probe_rotation(low_dec_cued_ind,1)==p);
        
        num_imp_uncued_high(sub,c)=sum(resp_high_uncued(probe_rotation(high_dec_uncued_ind,1)==p,1),1);
        num_imp_uncued_low(sub,c)=sum(resp_low_uncued(probe_rotation(low_dec_uncued_ind,1)==p,1),1);
        tot_imp_uncued_high(sub,c)=sum(probe_rotation(high_dec_uncued_ind,1)==p);
        tot_imp_uncued_low(sub,c)=sum(probe_rotation(low_dec_uncued_ind,1)==p);
        c=c+1;
    end
    probe_rot_levels=unique(probe_rotation)';
end
%% cluster corrected stats
[datobs, datrnd] = cluster_test_helper(dec_imp_cued(:,time>0)', 50000);
[h_imp_cued, p_imp_cued, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_imp_uncued(:,time>0)', 50000);
[h_imp_uncued, p_imp_uncued, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_imp_cued(:,time>0)'-dec_imp_uncued(:,time>0)', 50000);
[h_imp_diff, p_imp_diff, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

ci_imp_cued=bootci(50000,@mean,(dec_imp_cued)); % get C.I. for plotting
ci_imp_uncued=bootci(50000,@mean,(dec_imp_uncued));
%% plot decoding time-course and significant time-clusters

test_time=time(time>0);
pclustu_cued = unique(p_imp_cued);
npclust_cued = nnz(pclustu_cued < 0.05);
pclustu_uncued = unique(p_imp_uncued);
npclust_uncued = nnz(pclustu_uncued < 0.05);
pclustu_diff = unique(p_imp_diff);
npclust_diff = nnz(pclustu_diff < 0.05);

figure(1)
title('Figure 2F')
hold all
plot(time,mean(dec_imp_uncued,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp_cued,1),'Color',[0 0 1 1],'LineWidth',2)
line('XData', [-0.1 .5], 'YData', [0 0], 'LineStyle', '-','LineWidth', 1, 'Color','k');
fill([0,0,.1,.1],[-0.0005,-0.00035,-0.00035,-0.0005],[0 0 0],'EdgeColor','none')
plot(time,mean(dec_imp_uncued,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp_cued,1),'Color',[0 0 1 1],'LineWidth',2)
plot(time,ci_imp_uncued,'Color',[1 0 0 1],'LineWidth',.5, 'LineStyle', '-.')
plot(time,ci_imp_cued,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', '-.')
% extract time range of each significant cluster of each test and show in figure
for ipclust = 1:npclust_uncued
    currind  = p_imp_uncued == pclustu_uncued(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0022,0.0023,0.0023,0.0022],[1 0 0],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0023,0.0023,0],[1 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_cued
    currind  = p_imp_cued == pclustu_cued(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0024,0.0025,0.0025,0.0024],[0 0 1],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0025,0.0025,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_diff
    currind  = p_imp_diff == pclustu_diff(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0023,0.0024,0.0024,0.0023],[0.5 0 0.5],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0024,0.0024,0],[0.5 0 0.5],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-0.1 .5]);ylim([-0.0005 0.0025]);
set(gca,'TickDir','out')
xticks([linspace(-.1,.9,10)])
xlabel('Time (s) - relative to impulse')
ylabel('Decoding accuracy')
legend('uncued','cued','Location','northwest')
pbaspect([.85,1,1])

%% stats on averaged decoding
dec_imp_cued_av=mean(dec_imp_cued(:,time>.1),2);
dec_imp_uncued_av=mean(dec_imp_uncued(:,time>.1),2);

p_cued_av=GroupPermTest(dec_imp_cued_av,50000,2);
p_uncued_av=GroupPermTest(dec_imp_uncued_av,50000,2);
p_diff_av=GroupPermTest(dec_imp_cued_av-dec_imp_uncued_av,50000,2);

ci_imp_cued_av=bootci(50000,@mean,dec_imp_cued_av);
ci_imp_uncued_av=bootci(50000,@mean,dec_imp_uncued_av);
%% plot boxplots and error bars
figure(2)
title('Figure 2F boxplots')
hold all
b=boxplot([dec_imp_cued_av,dec_imp_uncued_av],[1 2]);
set(b(:,1),'color','b');
set(b(:,2),'color','r');
plot(1,mean(dec_imp_cued_av,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(2,mean(dec_imp_uncued_av,1),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_imp_cued_av','Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp_cued_av(1) ci_imp_cued_av(1) ],'Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp_cued_av(2) ci_imp_cued_av(2) ],'Color','b','LineWidth',3)
plot([2 2],ci_imp_uncued_av','Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp_uncued_av(1) ci_imp_uncued_av(1) ],'Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp_uncued_av(2) ci_imp_uncued_av(2) ],'Color','r','LineWidth',3)
ylim([-.0012 .0025])
%% stats on difference between high and low decoding trials
p_acc_cued=GroupPermTest(acc_diff_high_low_cued,50000,2);
p_acc_uncued=GroupPermTest(acc_diff_high_low_uncued,50000,2);
%% plotting accuracy difference between high and low decoding trials
ci_acc_cued=bootci(50000,@mean,acc_diff_high_low_cued);
ci_acc_uncued=bootci(50000,@mean,acc_diff_high_low_uncued);

%%
figure(3)
title('Figure 3A left')
boxplot(acc_diff_high_low_cued,'Color',[0 0 1]);
hold all
plot(1,mean(acc_diff_high_low_cued,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_acc_cued','Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_cued(1) ci_acc_cued(1) ],'Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_cued(2) ci_acc_cued(2) ],'Color',[0 0 1],'LineWidth',3)
ylim([-.08 .08])
ylabel('Mean accuracy difference')
title('cued item')
set(gca,'TickDir','out')

figure(4)
title('Figure 3B left')
boxplot(acc_diff_high_low_uncued,'Color',[1 0 0]);
hold all
plot(1,mean(acc_diff_high_low_uncued,1),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_acc_uncued','Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_uncued(1) ci_acc_uncued(1) ],'Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_uncued(2) ci_acc_uncued(2) ],'Color',[1 0 0],'LineWidth',3)
ylim([-.08 .08])
ylabel('Mean accuracy difference')
title('uncued item')
set(gca,'TickDir','out')
%% beahioural modelling...
StimLevels=[probe_rot_levels;probe_rot_levels];
params = [0 .2 0 0];
PF = @PAL_CumulativeNormal;

options = PAL_minimize('options');   %PAL_minimize search options
options.TolX = 1e-15;       %Increase desired precision on parameters
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 1];
% ... for high and low cued deocding trials
for sub=1:30
    fprintf(['Doing ' num2str(sub) '\n'])
    NumPos=[num_imp_cued_high(sub,:);num_imp_cued_low(sub,:)];
    OutOfNum=[tot_imp_cued_high(sub,:);tot_imp_cued_low(sub,:)];
    parameters = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
        OutOfNum, params, PF, 'thresholds','constrained','slopes',...
        'unconstrained','gammaEQlambda',1,'lapserates','unconstrained','guessrates','unconstrained',...
        'lapseLimits', lapseLimits,'SearchOptions',options);
    slope_cued_high(sub,1)=parameters(1,2);
    slope_cued_low(sub,1)=parameters(2,2);
    as_cued_high(sub,1)=parameters(1,3);
    as_cued_low(sub,1)=parameters(2,3);
    params_cued_high(sub,:)=parameters(1,:);
    params_cued_low(sub,:)=parameters(2,:);
end
p_slope_cued=GroupPermTest(slope_cued_high-slope_cued_low,50000,1); % one-sided test on the slope/precision difference
p_as_cued=GroupPermTest(as_cued_high-as_cued_low,50000,-1); % one-sided test on the asymptote/guess-rate difference
%% ... for high and low uncued deocding trials
for sub=1:30
    fprintf(['Doing ' num2str(sub) '\n'])
    NumPos=[num_imp_uncued_high(sub,:);num_imp_uncued_low(sub,:)];
    OutOfNum=[tot_imp_uncued_high(sub,:);tot_imp_uncued_low(sub,:)];
    parameters = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
        OutOfNum, params, PF, 'thresholds','constrained','slopes',...
        'unconstrained','gammaEQlambda',1,'lapserates','unconstrained','guessrates','unconstrained',...
        'lapseLimits', lapseLimits,'SearchOptions',options);
    slope_uncued_high(sub,1)=parameters(1,2);
    slope_uncued_low(sub,1)=parameters(2,2);
    as_uncued_high(sub,1)=parameters(1,3);
    as_uncued_low(sub,1)=parameters(2,3);
    params_uncued_high(sub,:)=parameters(1,:);
    params_uncued_low(sub,:)=parameters(2,:);
end
p_slope_uncued=GroupPermTest(slope_uncued_high-slope_uncued_low,50000,-1); % one-sided test on the slope/precision difference
p_as_uncued=GroupPermTest(as_uncued_high-as_uncued_low,50000,1); % one-sided test on the asymptote/guess-rate difference
%% plot sigmoids of models and data of cued item
StimLevelsFineGrain = [min(min(StimLevels)):(max(max(StimLevels) - ... 
    min(min(StimLevels))))./1000:max(max(StimLevels))];

for sub=1:30
    Model_cued_high(sub,:) = PF(params_cued_high(sub,:),StimLevelsFineGrain);
    Model_cued_low(sub,:) = PF(params_cued_low(sub,:),StimLevelsFineGrain);
    Model_uncued_high(sub,:) = PF(params_uncued_high(sub,:),StimLevelsFineGrain);
    Model_uncued_low(sub,:) = PF(params_uncued_low(sub,:),StimLevelsFineGrain);
end

ci_resp_cued_high=bootci(50000,@mean,num_imp_cued_high./tot_imp_cued_high);
ci_resp_cued_low=bootci(50000,@mean,num_imp_cued_low./tot_imp_cued_low);

ci_resp_uncued_high=bootci(50000,@mean,num_imp_uncued_high./tot_imp_uncued_high);
ci_resp_uncued_low=bootci(50000,@mean,num_imp_uncued_low./tot_imp_uncued_low);


figure(5)
title('Figure 3A right')
hold all
h=plot(probe_rot_levels,mean(num_imp_cued_high./tot_imp_cued_high,1),'o','Color',[0 0 .5]);
 set(h,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 .5])
plot(probe_rot_levels,mean(num_imp_cued_low./tot_imp_cued_low,1),'o','Color',[.4 .4 .8])
plot(StimLevelsFineGrain,mean(Model_cued_high,1),'-','linewidth',2,...
    'color',[0 0 .5]);
plot(StimLevelsFineGrain,mean(Model_cued_low,1),'-','linewidth',2,...
    'color',[.4 .4 .8]);
% plot([probe_rot_levels;probe_rot_levels],ci_resp_cued_high,'.','Color',[0 .5 0])
for p=1:length(probe_rot_levels)
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_cued_low(1,p) ci_resp_cued_low(1,p)],'Color',[0 0 .5])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_cued_low(2,p) ci_resp_cued_low(2,p)],'Color',[0 0 .5])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_cued_low(:,p)','Color',[0 0 .5])
    
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_cued_high(1,p) ci_resp_cued_high(1,p)],'Color',[.4 .4 .8])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_cued_high(2,p) ci_resp_cued_high(2,p)],'Color',[.4 .4 .8])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_cued_high(:,p)','Color',[.4 .4 .8])
end
legend('data high','data low','model high','model low','Location','northwest')
xlabel('degree rotation')
ylabel('proportion clockwise response')
title('cued item')
set(gca,'TickDir','out');ax = gca;ax.XTick = probe_rot_levels;
xlim([-43 43])
ylim([0 1])

ci_slope_cued=bootci(50000,@mean,slope_cued_high-slope_cued_low);
%% plot boxplot inset of cued precision difference
figure(6)
title('Figure 3A right inset')
hold all
boxplot(slope_cued_high-slope_cued_low,'Color',[0 0 1]);
plot(1,mean(slope_cued_high-slope_cued_low,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_slope_cued','Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_cued(1) ci_slope_cued(1) ],'Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_cued(2) ci_slope_cued(2) ],'Color',[0 0 1],'LineWidth',3)
ylim([-.08 .08])
ylabel('Precision difference')
title('cued item')
set(gca,'TickDir','out')

%% plot sigmoids of models and data of uncued item
figure(7)
title('Figure 3B right')
hold all
h=plot(probe_rot_levels,mean(num_imp_uncued_high./tot_imp_uncued_high,1),'o','Color',[.5 0 0]);
 set(h,'MarkerEdgeColor','none','MarkerFaceColor',[.5 0 0])
plot(probe_rot_levels,mean(num_imp_uncued_low./tot_imp_uncued_low,1),'o','Color',[.8 .4 .4])
plot(StimLevelsFineGrain,mean(Model_uncued_high,1),'-','linewidth',2,...
    'color',[.5 0 0]);
plot(StimLevelsFineGrain,mean(Model_uncued_low,1),'-','linewidth',2,...
    'color',[.8 .4 .4]);
for p=1:length(probe_rot_levels)
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_uncued_low(1,p) ci_resp_uncued_low(1,p)],'Color',[.5 0 0])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_uncued_low(2,p) ci_resp_uncued_low(2,p)],'Color',[.5 0 0])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_uncued_low(:,p)','Color',[.5 0 0])
    
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_uncued_high(1,p) ci_resp_uncued_high(1,p)],'Color',[.8 .4 .4])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_uncued_high(2,p) ci_resp_uncued_high(2,p)],'Color',[.8 .4 .4])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_uncued_high(:,p)','Color',[.8 .4 .4])
end
legend('data high','data low','model high','model low','Location','northwest')
xlabel('degree rotation')
ylabel('proportion clockwise response')
title('uncued item')
set(gca,'TickDir','out');ax = gca;ax.XTick = probe_rot_levels;
xlim([-43 43])
ylim([0 1])

%% plot boxplot inset of uncued precision difference
ci_slope_uncued=bootci(50000,@mean,slope_uncued_high-slope_uncued_low);

figure(8)
title('Figure 3B right inset')
hold all
boxplot(slope_uncued_high-slope_uncued_low,'Color',[0 0 1]);
plot(1,mean(slope_uncued_high-slope_uncued_low,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_slope_uncued','Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_uncued(1) ci_slope_uncued(1) ],'Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_uncued(2) ci_slope_uncued(2) ],'Color',[0 0 1],'LineWidth',3)
ylim([-.08 .08])
ylabel('Precision difference')
title('uncued item')
set(gca,'TickDir','out')


