% Script containing analyses for Figure 2d of "Dynamic hidden
% states underlying working-memory-guided behavior", Nature Neuroscience,
% 2017

% Input is preprocessed (see article) and baselined (-200 to 0 ms relative
% to stimulus onset) EEG data.

close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states'; % path to main folder containing all data, functions, toolboxes, scripts...
addpath(genpath(main_dir))
dat_dir=[main_dir '/Data'];

angspace=(-pi:pi/6:pi)'; % angular space for tuning curve (in radians)
angspace(end)=[];
bin_width=pi/6; % width of each angle bin of the tuning curve (in radians)

s_factor=8; % determines the SD of the smoothing kernel for the decoding time-courses,
            % depends on resolution of data (500 hz = 2ms, 2ms times 8 is thus an SD of 16 ms)

for sub=1:30
    fprintf(['Doing ' num2str(sub) '\n'])
    
    load(fullfile(dat_dir,['Dynamic_hidden_states_exp1_' num2str(sub) '.mat']));
    
    EEG_dat=exp1_data.EEG_mem_items;
    Results=exp1_data.Results;
    time=EEG_dat.time;
    
    clear exp1_data
    
    incl=not(ismember(1:size(EEG_dat.trial,1),EEG_dat.bad_trials))'; % logical array of trials to be included
    
    data=EEG_dat.trial(incl,:,:);
    
    clear EEG_dat
        
    data= bsxfun(@minus, data, mean(data,2)); % mean center voltage across channels to normalize
    
    mem_angles=Results(incl,1:2)*2; % extract memory item angles and rescale 
     
    dec_left = mahalTune_func(data,mem_angles(:,1),angspace,bin_width); % decode item presented on the left
    dec_right = mahalTune_func(data,mem_angles(:,2),angspace,bin_width); % decode item presented on the right
    
    dec_mem(sub,:)=imgaussfilt(mean(cat(1,dec_left,dec_right),1),s_factor); % average over decoding values of each item and all trials and smooth over time
end
%% cluster corrected statistics
[datobs, datrnd] = cluster_test_helper(dec_mem(:,time>0)', 50000);
[h_mem, p_mem, ~] = cluster_test(datobs,datrnd,0,0.05,0.05);

ci_mem=bootci(50000,@mean,(dec_mem)); % compute C.I. for plotting
%%
test_time=time(time>0);
pclustu = unique(p_mem);
npclust = nnz(pclustu < 0.05);
figure (1)
title('Figure 2D')
hold all
plot(time,mean(dec_mem,1),'Color',[.1 .1 .1 1],'LineWidth',2)
line('XData', [-0.1 1.05], 'YData', [0 0], 'LineStyle', '-','LineWidth', 1, 'Color','k');
fill([0,0,.25,.25],[-0.0005,-0.00035,-0.00035,-0.0005],[0.5 0.5 0.5],'EdgeColor','none')
plot(time,mean(dec_mem,1),'Color',[.1 .1 .1 1],'LineWidth',2)
plot(time,ci_mem(:,:),'Color',[.1 .1 .1 1],'LineWidth',.5, 'LineStyle', '-.')
for ipclust = 1:npclust % extract time range of each significant cluster and show in figure
    currind  = p_mem == pclustu(ipclust);
    fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0.0035,0.0036,0.0036,0.0035],[0 0 0.8],'EdgeColor','none')
    h=fill([min(test_time(currind)),min(test_time(currind)),max(test_time(currind)),max(test_time(currind))],[0,0.0036,0.0036,0],[0 0 0.8],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-0.1 1.05]);ylim([-0.0005 0.0036])
set(gca,'TickDir','out')
xlabel('Time (s) - relative to memory items')
legend('both items')