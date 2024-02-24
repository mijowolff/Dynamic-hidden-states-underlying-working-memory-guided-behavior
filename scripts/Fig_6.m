% Script containing analyses for Figures 6 of "Dynamic hidden
% states underlying working-memory-guided behavior", Nature Neuroscience,
% 2017

% Input is preprocessed (see article) and baselined (-200 to 0 ms relative
% to stimulus onset) EEG data.

close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states'; % path to main folder containing all data, functions/toolboxes, scripts...
addpath(genpath(main_dir))
dat_dir=  '/Data/Michael/Dynamic_hidden_states';

angspace=(-pi:pi/6:pi)'; % angular space for tuning curve (in radians)
angspace(end)=[];
bin_width=pi/6; % width of each angle bin of the tuning curve (in radians)

s_factor=8; % determines the SD of the smoothing kernel for the decoding time-courses,
% depends on resolution of data (500 hz = 2ms, 2ms times 8 is thus an SD of 16 ms)
%%
for sub=1:19
    fprintf(['Doing ' num2str(sub) '\n'])
    
    load(fullfile(dat_dir,['Dynamic_hidden_states_exp2_' num2str(sub) '.mat']));
    
    EEG_imp1_sess1=exp2_data.EEG_impulse1_sess1;
    EEG_imp2_sess1=exp2_data.EEG_impulse2_sess1;
    Results1=exp2_data.Results_sess1;
    
    EEG_imp1_sess2=exp2_data.EEG_impulse1_sess2;
    EEG_imp2_sess2=exp2_data.EEG_impulse2_sess2;
    Results2=exp2_data.Results_sess2;
    
    time=EEG_imp1_sess1.time;
    
    clear exp2_data
    
    incl_imp1_1=not(ismember(1:size(EEG_imp1_sess1.trial,1),EEG_imp1_sess1.bad_trials))'; % trials to be included
    incl_imp2_1=not(ismember(1:size(EEG_imp2_sess1.trial,1),EEG_imp2_sess1.bad_trials))';
    data_imp1_sess1 = EEG_imp1_sess1.trial(incl_imp1_1,:,:);
    data_imp2_sess1 = EEG_imp2_sess1.trial(incl_imp2_1,:,:);
    
    clear EEG_imp1_sess1 EEG_imp2_sess1
    
    incl_imp1_2=not(ismember(1:length(EEG_imp1_sess2.trial),EEG_imp1_sess2.bad_trials))'; % trials to be included
    incl_imp2_2=not(ismember(1:length(EEG_imp2_sess2.trial),EEG_imp2_sess2.bad_trials))';
    data_imp1_sess2 = EEG_imp1_sess2.trial(incl_imp1_2,:,:);
    data_imp2_sess2 = EEG_imp2_sess2.trial(incl_imp2_2,:,:);
    
    clear EEG_imp1_sess2 EEG_imp2_sess2
    
    data_imp1_sess1= bsxfun(@minus, data_imp1_sess1, mean(data_imp1_sess1,2)); % mean center voltage across channels (sort of normalization)
    data_imp2_sess1= bsxfun(@minus, data_imp2_sess1, mean(data_imp2_sess1,2));
    mem_angles1=Results1(:,1:2)*2; % extract memory item anglesand rescale
    
    data_imp1_sess2= bsxfun(@minus, data_imp1_sess2, mean(data_imp1_sess2,2)); % mean center voltage across channels (sort of normalization)
    data_imp2_sess2= bsxfun(@minus, data_imp2_sess2, mean(data_imp2_sess2,2));
    mem_angles2=Results2(:,1:2)*2; % extract memory item angles and rescale
    
    dec_imp1_early1 = mahalTune_func(data_imp1_sess1,mem_angles1(incl_imp1_1,1),angspace,bin_width); % decode early-tested item at impulse 1 in session 1
    dec_imp1_late1 = mahalTune_func(data_imp1_sess1,mem_angles1(incl_imp1_1,2),angspace,bin_width); % decode late-tested item at impulse 1 in session 1
    dec_imp2_early1 = mahalTune_func(data_imp2_sess1,mem_angles1(incl_imp2_1,1),angspace,bin_width); % decode early-tested item at impulse 2 in session 1
    dec_imp2_late1 = mahalTune_func(data_imp2_sess1,mem_angles1(incl_imp2_1,2),angspace,bin_width); % decode late-tested item at impulse 2 in session 1
    
    dec_imp1_early2 = mahalTune_func(data_imp1_sess2,mem_angles2(incl_imp1_2,1),angspace,bin_width); % decode early-tested item at impulse 1 in session 2
    dec_imp1_late2 = mahalTune_func(data_imp1_sess2,mem_angles2(incl_imp1_2,2),angspace,bin_width); % decode late-tested item at impulse 1 in session 2
    dec_imp2_early2 = mahalTune_func(data_imp2_sess2,mem_angles2(incl_imp2_2,1),angspace,bin_width); % decode early-tested item at impulse 2 in session 2
    dec_imp2_late2 = mahalTune_func(data_imp2_sess2,mem_angles2(incl_imp2_2,2),angspace,bin_width); % decode late-tested item at impulse 2 in session 2
    
    dec_imp1_early(sub,:)=imgaussfilt((mean(dec_imp1_early1,1)+mean(dec_imp1_early2,1))/2,s_factor); % average over trials and sessions, and smooth over time
    dec_imp1_late(sub,:)=imgaussfilt((mean(dec_imp1_late1,1)+mean(dec_imp1_late2,1))/2,s_factor);
    
    dec_imp2_early(sub,:)=imgaussfilt((mean(dec_imp2_early1,1)+mean(dec_imp2_early2,1))/2,s_factor);
    dec_imp2_late(sub,:)=imgaussfilt((mean(dec_imp2_late1,1)+mean(dec_imp2_late2,1))/2,s_factor);
    
    %% the following section extracts high and low decoding trials for figure 6b
    
    early_acc1=Results1(incl_imp1_1,6);
    early_acc2=Results2(incl_imp1_2,6);
    late_acc1=Results1(incl_imp2_1,7);
    late_acc2=Results2(incl_imp2_2,7);
    
    early_pr_rotation=cat(1,Results1(incl_imp1_1,4),Results2(incl_imp1_2,4));
    late_pr_rotation=cat(1,Results1(incl_imp2_1,5),Results2(incl_imp2_2,5));
    
    early_pr_rotation([isnan(early_acc1);isnan(early_acc2)])=[]; %remove non-responses
    
    dec_imp1_early1_av=mean(dec_imp1_early1(~isnan(early_acc1),time>.1),2); %average over time and exclude all non-responses
    dec_imp2_late1_av=mean(dec_imp2_late1(:,time>.1),2);
    dec_imp1_early2_av=mean(dec_imp1_early2(~isnan(early_acc2),time>.1),2);
    dec_imp2_late2_av=mean(dec_imp2_late2(:,time>.1),2);
    
    early_acc1(isnan(early_acc1))=[];
    early_acc2(isnan(early_acc2))=[];
    
    high_imp1_early1=cat(1,dec_imp1_early1_av>median(dec_imp1_early1_av)); % median-split high and low decoding trials
    high_imp1_early2=cat(1,dec_imp1_early2_av>median(dec_imp1_early2_av));
    low_imp1_early1=cat(1,dec_imp1_early1_av<median(dec_imp1_early1_av)); 
    low_imp1_early2=cat(1,dec_imp1_early2_av<median(dec_imp1_early2_av));
    
    high_imp2_late1=cat(1,dec_imp2_late1_av>median(dec_imp2_late1_av));
    high_imp2_late2=cat(1,dec_imp2_late2_av>median(dec_imp2_late2_av));
    low_imp2_late1=cat(1,dec_imp2_late1_av<median(dec_imp2_late1_av));
    low_imp2_late2=cat(1,dec_imp2_late2_av<median(dec_imp2_late2_av));
    
    
    acc_early_diff_imp1(sub,1)=(mean(early_acc1(high_imp1_early1,1),1)+mean(early_acc2(high_imp1_early2,1),1))/2 ... % take difference between accuracy of high
        -(mean(early_acc1(low_imp1_early1,1),1)+mean(early_acc2(low_imp1_early2,1),1))/2;                             %...and low decoding trials
    
    acc_late_diff_imp2(sub,1)=(mean(late_acc1(high_imp2_late1,1),1)+mean(late_acc2(high_imp2_late2,1),1))/2 ...
        -(mean(late_acc1(low_imp2_late1,1),1)+mean(late_acc2(low_imp2_late2,1),1))/2;
 
    %% convert accuracy to cw (1) or acw (0) response for behavioural modelling later on

    early_resp = [early_acc1;early_acc2]; early_resp(early_resp==0)=-1;early_resp = early_resp.* sign(early_pr_rotation);early_resp(early_resp==-1)=0;
    late_resp = [late_acc1;late_acc2]; late_resp(late_resp==0)=-1;late_resp = late_resp.* sign(late_pr_rotation);late_resp(late_resp==-1)=0;
    
    early_resp_high=early_resp([high_imp1_early1;high_imp1_early2],1);
    early_resp_low=early_resp([low_imp1_early1;low_imp1_early2],1);
    late_resp_high=late_resp([high_imp2_late1;high_imp2_late2],1);
    late_resp_low=late_resp([low_imp2_late1;low_imp2_late2],1);
    
    % prepare for modelling
    c=1;
    for p=unique(early_pr_rotation)'
        num_early_resp_high(sub,c)=sum(early_resp_high(early_pr_rotation([high_imp1_early1;high_imp1_early2],1)==p,1),1); % number of clockwise responses for each probe-rotation condition separately
        num_early_resp_low(sub,c)=sum(early_resp_low(early_pr_rotation([low_imp1_early1;low_imp1_early2],1)==p,1),1);
        tot_early_high(sub,c)=sum(early_pr_rotation([high_imp1_early1;high_imp1_early2],1)==p); % total number of specific probe-rotaion condition
        tot_early_low(sub,c)=sum(early_pr_rotation([low_imp1_early1;low_imp1_early2],1)==p);
        
        num_late_resp_high(sub,c)=sum(late_resp_high(late_pr_rotation([high_imp2_late1;high_imp2_late2],1)==p,1),1); % number of clockwise responses for each probe-rotation condition separately
        num_late_resp_low(sub,c)=sum(late_resp_low(late_pr_rotation([low_imp2_late1;low_imp2_late2],1)==p,1),1);
        tot_late_high(sub,c)=sum(late_pr_rotation([high_imp2_late1;high_imp2_late2],1)==p); % total number of specific probe-rotaion condition
        tot_late_low(sub,c)=sum(late_pr_rotation([low_imp2_late1;low_imp2_late2],1)==p);
        c=c+1;
    end
    probe_rot_levels=unique(early_pr_rotation)';
end
%% cluster corrected stats
[datobs, datrnd] = cluster_test_helper(dec_imp1_early(:,time>0)', 50000);
[h_imp1_early, p_imp1_early, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_imp1_late(:,time>0)', 50000);
[h_imp1_late, p_imp1_late, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

ci_imp1_early=bootci(50000,@mean,(dec_imp1_early)); % get C.I. for plotting
ci_imp1_late=bootci(50000,@mean,(dec_imp1_late));

[datobs, datrnd] = cluster_test_helper(dec_imp2_early(:,time>0)', 50000);
[h_imp2_early, p_imp2_early, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_imp2_late(:,time>0)', 50000);
[h_imp2_late, p_imp2_late, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

ci_imp2_early=bootci(50000,@mean,(dec_imp2_early)); % get C.I. for plotting
ci_imp2_late=bootci(50000,@mean,(dec_imp2_late));
%% plot decoding time-course and significant time-clusters
test_time=time(time>0);
pclustu_imp1_early = unique(p_imp1_early);
npclust_imp1_early = nnz(pclustu_imp1_early < 0.05);
pclustu_imp1_late = unique(p_imp1_late);
npclust_imp1_late = nnz(pclustu_imp1_late < 0.05);

figure(1)
title('Figure 6A') 
hold all
plot(time,mean(dec_imp1_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp1_early,1),'Color',[0 0 1 1],'LineWidth',2)
line('XData', [-0.1 .5], 'YData', [0 0], 'LineStyle', '-','LineWidth', 1, 'Color','k');
fill([0,0,.1,.1],[-0.0005,-0.00035,-0.00035,-0.0005],[0 0 0],'EdgeColor','none')
plot(time,mean(dec_imp1_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp1_early,1),'Color',[0 0 1 1],'LineWidth',2)
plot(time,ci_imp1_late,'Color',[1 0 0 1],'LineWidth',.5, 'LineStyle', '-.')
plot(time,ci_imp1_early,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', '-.')
% extract time range of each significant cluster of each test and show in figure
for ipclust = 1:npclust_imp1_late
    currind  = p_imp1_late == pclustu_imp1_late(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0022,0.0023,0.0023,0.0022],[1 0 0],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0023,0.0023,0],[1 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_imp1_early
    currind  = p_imp1_early == pclustu_imp1_early(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0024,0.0025,0.0025,0.0024],[0 0 1],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0025,0.0025,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-0.1 .5]);ylim([-0.0005 0.0025]);
set(gca,'TickDir','out')
xlabel('Time (s) - relative to impulse')
ylabel('Decoding accuracy')
legend('uncued','cued','Location','northwest')
pbaspect([.7,1,1])
%%
test_time=time(time>0);
pclustu_imp2_early = unique(p_imp2_early);
npclust_imp2_early = nnz(pclustu_imp2_early < 0.05);
pclustu_imp2_late = unique(p_imp2_late);
npclust_imp2_late = nnz(pclustu_imp2_late < 0.05);

figure(2)
title('Figure 6B') 
hold all
plot(time,mean(dec_imp2_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp2_early,1),'Color',[0 0 1 1],'LineWidth',2)
line('XData', [-0.1 .5], 'YData', [0 0], 'LineStyle', '-','LineWidth', 1, 'Color','k');
fill([0,0,.1,.1],[-0.0005,-0.00035,-0.00035,-0.0005],[0 0 0],'EdgeColor','none')
plot(time,mean(dec_imp2_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_imp2_early,1),'Color',[0 0 1 1],'LineWidth',2)
plot(time,ci_imp2_late,'Color',[1 0 0 1],'LineWidth',.5, 'LineStyle', '-.')
plot(time,ci_imp2_early,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', '-.')
% extract time range of each significant cluster of each test and show in figure
for ipclust = 1:npclust_imp2_late
    currind  = p_imp2_late == pclustu_imp2_late(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0022,0.0023,0.0023,0.0022],[1 0 0],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0023,0.0023,0],[1 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_imp2_early
    currind  = p_imp2_early == pclustu_imp2_early(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0024,0.0025,0.0025,0.0024],[0 0 1],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0025,0.0025,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-0.1 .5]);ylim([-0.0005 0.0025]);
set(gca,'TickDir','out')
xlabel('Time (s) - relative to impulse')
ylabel('Decoding accuracy')
legend('uncued','cued','Location','northwest')
pbaspect([.7,1,1])
%% stats on average decoding

dec_imp1_early_av=mean(dec_imp1_early(:,time>.1),2);
dec_imp1_late_av=mean(dec_imp1_late(:,time>.1),2);
dec_imp2_early_av=mean(dec_imp2_early(:,time>.1),2);
dec_imp2_late_av=mean(dec_imp2_late(:,time>.1),2);

p_imp1_early_av=GroupPermTest(dec_imp1_early_av,50000,2);
p_imp1_late_av=GroupPermTest(dec_imp1_late_av,50000,2);
p_imp1_diff_av=GroupPermTest(dec_imp1_early_av-dec_imp1_late_av,50000,2);

p_imp2_early_av=GroupPermTest(dec_imp2_early_av,50000,2);
p_imp2_late_av=GroupPermTest(dec_imp2_late_av,50000,2);
p_imp2_diff_av=GroupPermTest(dec_imp2_early_av-dec_imp2_late_av,50000,2);

ci_imp1_early_av=bootci(50000,@mean,dec_imp1_early_av);
ci_imp1_late_av=bootci(50000,@mean,dec_imp1_late_av);
ci_imp2_early_av=bootci(50000,@mean,dec_imp2_early_av);
ci_imp2_late_av=bootci(50000,@mean,dec_imp2_late_av);

%% plot boxplots of average decoding
figure(3)
title('Figure 6A, boxplots')
hold all
b=boxplot([dec_imp1_early_av,dec_imp1_late_av],[1 2]);
set(b(:,1),'color','b');
set(b(:,2),'color','r');
plot(1,mean(dec_imp1_early_av,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(2,mean(dec_imp1_late_av,1),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_imp1_early_av','Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp1_early_av(1) ci_imp1_early_av(1) ],'Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp1_early_av(2) ci_imp1_early_av(2) ],'Color','b','LineWidth',3)
plot([2 2],ci_imp1_late_av','Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp1_late_av(1) ci_imp1_late_av(1) ],'Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp1_late_av(2) ci_imp1_late_av(2) ],'Color','r','LineWidth',3)
ylim([-.0006 .0025])

figure(4)
title('Figure 6B, boxplots')
hold all
b=boxplot([dec_imp2_early_av,dec_imp2_late_av],[1 2]);
set(b(:,1),'color','b');
set(b(:,2),'color','r');
plot(1,mean(dec_imp2_early_av,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(2,mean(dec_imp2_late_av,1),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_imp2_early_av','Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp2_early_av(1) ci_imp2_early_av(1) ],'Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_imp2_early_av(2) ci_imp2_early_av(2) ],'Color','b','LineWidth',3)
plot([2 2],ci_imp2_late_av','Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp2_late_av(1) ci_imp2_late_av(1) ],'Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_imp2_late_av(2) ci_imp2_late_av(2) ],'Color','r','LineWidth',3)
ylim([-.0006 .0025])

%% plotting accuracy difference between high and low decoding trials
p_acc_imp1_early=GroupPermTest(acc_early_diff_imp1,50000,2);
p_acc_imp2_late=GroupPermTest(acc_late_diff_imp2,50000,2);
%% beahioural modelling
StimLevels=[probe_rot_levels;probe_rot_levels];
params = [0 .2 0 0];
PF = @PAL_CumulativeNormal;

options = PAL_minimize('options');   %PAL_minimize search options
options.TolX = 1e-15;       %Increase desired precision on parameters
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 1];

for sub=1:19
    fprintf(['Doing ' num2str(sub) '\n'])
    NumPos=[num_early_resp_high(sub,:);num_early_resp_low(sub,:)];
    OutOfNum=[tot_early_high(sub,:);tot_early_low(sub,:)];
    parameters = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
        OutOfNum, params, PF, 'thresholds','constrained','slopes',...
        'unconstrained','gammaEQlambda',1,'lapserates','unconstrained','guessrates','unconstrained',...
        'lapseLimits', lapseLimits,'SearchOptions',options);
    slope_early_high(sub,1)=parameters(1,2);
    slope_early_low(sub,1)=parameters(2,2);
    as_early_high(sub,1)=parameters(1,3);
    as_early_low(sub,1)=parameters(2,3);
    params_early_high(sub,:)=parameters(1,:);
    params_early_low(sub,:)=parameters(2,:);
end

p_slope_early=GroupPermTest(slope_early_high-slope_early_low,50000,1); % one-sided test on the slope/precision difference
p_as_early=GroupPermTest(as_early_high-as_early_low,50000,-1);% one-sided test on the asymptote/guess-rate difference

for sub=1:19
    fprintf(['Doing ' num2str(sub) '\n'])
    NumPos=[num_late_resp_high(sub,:);num_late_resp_low(sub,:)];
    OutOfNum=[tot_late_high(sub,:);tot_late_low(sub,:)];
    parameters = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
        OutOfNum, params, PF, 'thresholds','constrained','slopes',...
        'unconstrained','gammaEQlambda',1,'lapserates','unconstrained','guessrates','unconstrained',...
        'lapseLimits', lapseLimits,'SearchOptions',options);
    slope_late_high(sub,1)=parameters(1,2);
    slope_late_low(sub,1)=parameters(2,2);
    as_late_high(sub,1)=parameters(1,3);
    as_late_low(sub,1)=parameters(2,3);
    params_late_high(sub,:)=parameters(1,:);
    params_late_low(sub,:)=parameters(2,:);
end

p_slope_late=GroupPermTest(slope_late_high-slope_late_low,50000,1); % one-sided test on the slope/precision difference
p_as_late=GroupPermTest(as_late_high-as_late_low,50000,-1); % one-sided test on the asymptote/guess-rate difference
%% plot accuracy and model differences
ci_acc_early=bootci(50000,@mean,acc_early_diff_imp1);
ci_acc_late=bootci(50000,@mean,acc_late_diff_imp2);

figure(5)
title('Figure 6C left')
hold all
boxplot(acc_early_diff_imp1,'Color',[0 0 1]);
plot(1,mean(acc_early_diff_imp1,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_acc_early','Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_early(1) ci_acc_early(1) ],'Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_early(2) ci_acc_early(2) ],'Color',[0 0 1],'LineWidth',3)
ylim([-.08 .08])
ylabel('Mean accuracy difference')
title('early item')
set(gca,'TickDir','out')

figure(6)
title('Figure 6D left')
hold all
boxplot(acc_late_diff_imp2,'Color',[1 0 0]);
plot(1,mean(acc_late_diff_imp2,1),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_acc_late','Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_late(1) ci_acc_late(1) ],'Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_acc_late(2) ci_acc_late(2) ],'Color',[1 0 0],'LineWidth',3)
ylim([-.08 .08])
ylabel('Mean accuracy difference')
title('late item')
set(gca,'TickDir','out')

%%
StimLevelsFineGrain = [min(min(StimLevels)):(max(max(StimLevels) - ... 
    min(min(StimLevels))))./1000:max(max(StimLevels))];

for sub=1:19
    Model_early_high(sub,:) = PF(params_early_high(sub,:),StimLevelsFineGrain);
    Model_early_low(sub,:) = PF(params_early_low(sub,:),StimLevelsFineGrain);
    Model_late_high(sub,:) = PF(params_late_high(sub,:),StimLevelsFineGrain);
    Model_late_low(sub,:) = PF(params_late_low(sub,:),StimLevelsFineGrain);
end

ci_resp_early_high=bootci(50000,@mean,num_early_resp_high./tot_early_high);
ci_resp_early_low=bootci(50000,@mean,num_early_resp_low./tot_early_low);

ci_resp_late_high=bootci(50000,@mean,num_late_resp_high./tot_late_high);
ci_resp_late_low=bootci(50000,@mean,num_late_resp_low./tot_late_low);
%%
figure(7)
title('Figure 6C right')
hold all
h=plot(probe_rot_levels,mean(num_early_resp_high./tot_early_high,1),'o','Color',[0 0 .5]);
 set(h,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 .5])
plot(probe_rot_levels,mean(num_early_resp_low./tot_early_low,1),'o','Color',[.4 .4 .8])
plot(StimLevelsFineGrain,mean(Model_early_high,1),'-','linewidth',2,...
    'color',[0 0 .5]);
plot(StimLevelsFineGrain,mean(Model_early_low,1),'-','linewidth',2,...
    'color',[.4 .4 .8]);
for p=1:length(probe_rot_levels)
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_early_low(1,p) ci_resp_early_low(1,p)],'Color',[.4 .4 .8])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_early_low(2,p) ci_resp_early_low(2,p)],'Color',[.4 .4 .8])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_early_low(:,p)','Color',[.4 .4 .8])
    
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_early_high(1,p) ci_resp_early_high(1,p)],'Color',[0 0 .5])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_early_high(2,p) ci_resp_early_high(2,p)],'Color',[0 0 .5])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_early_high(:,p)','Color',[0 0 .5])
end
legend('data high','data low','model high','model low','Location','northwest')
xlabel('degree rotation')
ylabel('proportion clockwise response')
title('early item')
set(gca,'TickDir','out');ax = gca;ax.XTick = probe_rot_levels;
xlim([-43 43])
ylim([0 1])
%%
figure(8)
title('Figure 6D right')
hold all
h=plot(probe_rot_levels,mean(num_late_resp_high./tot_late_high,1),'o','Color',[.5 0 0]);
 set(h,'MarkerEdgeColor','none','MarkerFaceColor',[.5 0 0])
plot(probe_rot_levels,mean(num_late_resp_low./tot_late_low,1),'o','Color',[.8 .4 .4])
plot(StimLevelsFineGrain,mean(Model_late_high,1),'-','linewidth',2,...
    'color',[.5 0 0]);
plot(StimLevelsFineGrain,mean(Model_late_low,1),'-','linewidth',2,...
    'color',[.8 .4 .4]);
for p=1:length(probe_rot_levels)
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_late_low(1,p) ci_resp_late_low(1,p)],'Color',[.8 .4 .4])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_late_low(2,p) ci_resp_late_low(2,p)],'Color',[.8 .4 .4])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_late_low(:,p)','Color',[.8 .4 .4])
    
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_late_high(1,p) ci_resp_late_high(1,p)],'Color',[.5 0 0])
    plot([probe_rot_levels(p)-.5 probe_rot_levels(p)+.5],[ci_resp_late_high(2,p) ci_resp_late_high(2,p)],'Color',[.5 0 0])
    plot([probe_rot_levels(p) probe_rot_levels(p)],ci_resp_late_high(:,p)','Color',[.5 0 0])
end
legend('data high','data low','model high','model low','Location','northwest')
xlabel('degree rotation')
ylabel('proportion clockwise response')
title('late item')
set(gca,'TickDir','out');ax = gca;ax.XTick = probe_rot_levels;
xlim([-43 43])
ylim([0 1])
%%
ci_slope_early=bootci(50000,@mean,slope_early_high-slope_early_low);
ci_slope_late=bootci(50000,@mean,slope_late_high-slope_late_low);

figure(9)
title('Figure 6C right inset')
hold all
boxplot(slope_early_high-slope_early_low,'Color',[0 0 1]);
plot(1,mean(slope_early_high-slope_early_low,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_slope_early','Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_early(1) ci_slope_early(1) ],'Color',[0 0 1],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_early(2) ci_slope_early(2) ],'Color',[0 0 1],'LineWidth',3)
ylim([-.08 .08])
ylabel('Precision difference')
title('early item')
set(gca,'TickDir','out')

figure(10)
title('Figure 6D right inset')
hold all
boxplot(slope_late_high-slope_late_low,'Color',[1 0 0]);
plot(1,mean(slope_late_high-slope_late_low,1),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_slope_late','Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_late(1) ci_slope_late(1) ],'Color',[1 0 0],'LineWidth',3)
plot([1-.02 1+.02],[ci_slope_late(2) ci_slope_late(2) ],'Color',[1 0 0],'LineWidth',3)
ylim([-.04 .15])
ylabel('Precision difference')
title('late item')
set(gca,'TickDir','out')

