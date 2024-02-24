% Script containing analyses for Figure 5a of "Dynamic hidden
% states underlying working-memory-guided behavior", Nature Neuroscience,
% 2017

% Input is preprocessed (see article) and baselined (-200 to 0 ms relative
% to stimulus onset) EEG data.

close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states'; % path to main folder containing all data, functions, toolboxes, scripts
addpath(genpath(main_dir))
dat_dir=[main_dir '\Data'];

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
    
    data1= bsxfun(@minus, data1, mean(data1,2)); % mean center voltage across channels to normalize
    mem_angles1=Results1(incl1,1:2)*2; % extract memory item angles and rescale 
    
    data2= bsxfun(@minus, data2, mean(data2,2)); % mean center voltage across channels to normalize
    mem_angles2=Results2(incl2,1:2)*2; % extract memory item angles and rescale

    dec_early1 = mahalTune_func(data1,mem_angles1(:,1),angspace,bin_width); % decode early-tested item in session 1
    dec_late1 = mahalTune_func(data1,mem_angles1(:,2),angspace,bin_width); % decode late-tested item in session 1
    
    dec_early2 = mahalTune_func(data2,mem_angles2(:,1),angspace,bin_width); % decode early-tested item in session 2
    dec_late2 = mahalTune_func(data2,mem_angles2(:,2),angspace,bin_width); % decode late-tested item in session 2 
    
    dec_mem_early(sub,:)=imgaussfilt((mean(dec_early1,1)+mean(dec_early2,1))/2,s_factor); % average over decoding values of each item and all trials and smooth the time-course
    dec_mem_late(sub,:)=imgaussfilt((mean(dec_late1,1)+mean(dec_late2,1))/2,s_factor);       
end
%% cluster corrected stats
[datobs, datrnd] = cluster_test_helper(dec_mem_early(:,time>0)', 50000);
[h_mem_early, p_mem_early, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_mem_late(:,time>0)', 50000);
[h_mem_late,p_mem_late,~] = cluster_test(datobs,datrnd,0,0.05,0.05);

[datobs, datrnd] = cluster_test_helper(dec_mem_early(:,time>0)'-dec_mem_late(:,time>0)', 50000);
[h_mem_diff, p_mem_diff, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

ci_mem_early=bootci(50000,@mean,(dec_mem_early)); % get C.I. for plotting
ci_mem_late=bootci(50000,@mean,(dec_mem_late)); 

%% plot decoding time-course and significant time-clusters

test_time=time(time>0);
pclustu_early = unique(p_mem_early);
npclust_early = nnz(pclustu_early < 0.05);
pclustu_late = unique(p_mem_late);
npclust_late = nnz(pclustu_late < 0.05);
pclustu_diff = unique(p_mem_diff);
npclust_diff = nnz(pclustu_diff < 0.05);

figure(1);
title('Figure 5A')
hold all
plot(time,mean(dec_mem_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_mem_early,1),'Color',[0 0 1 1],'LineWidth',2)
line('XData', [-0.1 1.2], 'YData', [0 0], 'LineStyle', '-','LineWidth', 1, 'Color','k');
fill([0,0,.25,.25],[-0.0005,-0.00035,-0.00035,-0.0005],[.5 .5 .5],'EdgeColor','none')
plot(time,mean(dec_mem_late,1),'Color',[1 0 0 1],'LineWidth',2)
plot(time,mean(dec_mem_early,1),'Color',[0 0 1 1],'LineWidth',2)
plot(time,ci_mem_late,'Color',[1 0 0 1],'LineWidth',.5, 'LineStyle', '-.')
plot(time,ci_mem_early,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', '-.')
% extract time range of each significant cluster of each test and show in figure
for ipclust = 1:npclust_late 
    currind  = p_mem_late == pclustu_late(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0027,0.0028,0.0028,0.0027],[1 0 0],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0028,0.0028,0],[1 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_early 
    currind  = p_mem_early == pclustu_early(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0029,0.003,0.003,0.0029],[0 0 1],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.003,0.003,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_diff 
    currind  = p_mem_diff == pclustu_diff(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0028,0.0029,0.0029,0.0028],[0.5 0 0.5],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0029,0.0029,0],[0.5 0 0.5],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-0.1 1.2]);ylim([-0.0005 0.0031]);
set(gca,'TickDir','out')
xlabel('Time (s) - relative to impulse')
ylabel('Decoding accuracy')
legend('uncued','cued','Location','northeast')

%% stats on averaged decoding
dec_mem_early_av=mean(dec_mem_early(:,time>.1),2);
dec_mem_late_av=mean(dec_mem_late(:,time>.1),2);

p_early_av=GroupPermTest(dec_mem_early_av,50000,2);
p_late_av=GroupPermTest(dec_mem_late_av,50000,2);
p_diff_av=GroupPermTest(dec_mem_early_av-dec_mem_late_av,50000,2);

ci_mem_early_av=bootci(50000,@mean,dec_mem_early_av);
ci_mem_late_av=bootci(50000,@mean,dec_mem_late_av);

%% plot boxplots and error bars
figure(2)
title('Figure 5A boxplots')
hold all
b=boxplot([dec_mem_early_av,dec_mem_late_av],[1 2]);
set(b(:,1),'color','b');
set(b(:,2),'color','r');
plot(1,mean(dec_mem_early_av,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(2,mean(dec_mem_late_av,1),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10)
plot([1 1],ci_mem_early_av','Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_mem_early_av(1) ci_mem_early_av(1) ],'Color','b','LineWidth',3)
plot([1-.1 1+.1],[ci_mem_early_av(2) ci_mem_early_av(2) ],'Color','b','LineWidth',3)
plot([2 2],ci_mem_late_av','Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_mem_late_av(1) ci_mem_late_av(1) ],'Color','r','LineWidth',3)
plot([2-.1 2+.1],[ci_mem_late_av(2) ci_mem_late_av(2) ],'Color','r','LineWidth',3)
ylim([-.0005 .0032])
