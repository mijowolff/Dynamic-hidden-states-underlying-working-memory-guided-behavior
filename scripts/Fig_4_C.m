close all
clear all

main_dir='/Data/Michael/Dropbox/Wolff/Wolff2016/Dynamic_hidden_states';
addpath(genpath(main_dir))
% dat_dir=[main_dir '/Data'];
dat_dir=  '/Data/Michael/Dynamic_hidden_states';
ft_path='/home/michael/Documents/MATLAB/fieldtrip';
addpath(ft_path)
ft_defaults

left_ch = {'P7';'P5';'P3';'P1';'PO7';'PO3';'O1'};
right_ch = {'P8';'P6';'P4';'P2';'PO8';'PO4';'O2'};
%%
for sub=1:19
    fprintf(['Doing ' num2str(sub) '\n'])
    
    load(fullfile(dat_dir,['Dynamic_hidden_states_exp2_whole_' num2str(sub) '.mat']));
        
    toi=find(exp2_data.EEG_mem_items_sess1.time>=-0.1&exp2_data.EEG_mem_items_sess1.time<=4.8);
    time=exp2_data.EEG_mem_items_sess1.time;
    
    incl1=setdiff(1:size(exp2_data.EEG_mem_items_sess1.trial,1),exp2_data.EEG_mem_items_sess1.bad_trials);
    incl2=setdiff(1:size(exp2_data.EEG_mem_items_sess2.trial,1),exp2_data.EEG_mem_items_sess2.bad_trials);

    Results1=exp2_data.Results_sess1(incl1,:);
    Results2=exp2_data.Results_sess2(incl2,:);

    frequencies=[6:0.5:16];
    cfg              = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi = frequencies;
    cfg.t_ftimwin = 5./cfg.foi;   % length of time window is 5 cycles
    cfg.toi = time(toi(1:5:end));   % time window "slides" from -.1 to 4.8 sec in steps of 0.01 sec
    cfg.keeptrials = 'yes';
    warning('off','all');
    
    cfg.trials      = incl1;
    freq_mem_1 = ft_freqanalysis(cfg, exp2_data.EEG_mem_items_sess1); % run for first
    cfg.trials      = incl2;
    freq_mem_2 = ft_freqanalysis(cfg, exp2_data.EEG_mem_items_sess2); % and second session
    clear exp2_data
    
    chan_labels=freq_mem_1.label;
    
    data1=log10(freq_mem_1.powspctrm);    
    clear freq_mem_1

    data2=log10(freq_mem_2.powspctrm);
    clear freq_mem_2
    
    if Results1(1,3)==1
        freq_lat_1=squeeze(mean(mean(data1(:,ismember(chan_labels,right_ch),:,:),1),2))-squeeze(mean(mean(data1(:,ismember(chan_labels,left_ch),:,:),1),2));
        freq_lat_2=squeeze(mean(mean(data2(:,ismember(chan_labels,left_ch),:,:),1),2))-squeeze(mean(mean(data2(:,ismember(chan_labels,right_ch),:,:),1),2));
        freq_lateralization(sub,:,:)=(freq_lat_1+freq_lat_2)/2;
    else 
        freq_lat_2=squeeze(mean(mean(data2(:,ismember(chan_labels,right_ch),:,:),1),2))-squeeze(mean(mean(data2(:,ismember(chan_labels,left_ch),:,:),1),2));
        freq_lat_1=squeeze(mean(mean(data1(:,ismember(chan_labels,left_ch),:,:),1),2))-squeeze(mean(mean(data1(:,ismember(chan_labels,right_ch),:,:),1),2));
        freq_lateralization(sub,:,:)=(freq_lat_1+freq_lat_2)/2;
    end
    clear data1 data2
end
%%
[datobs, datrnd] = cluster_test_helper(permute(freq_lateralization(:,:,:),[2 3 1]), 50000);
[h, pm_diff,~] = cluster_test(datobs,datrnd,0,0.05,0.05);
%%
B = bwboundaries(h);
fhandle=figure;
imagesc(time,frequencies,squeeze(mean(freq_lateralization,1)),[-0.05 0.05]); axis xy
hold all
for k=1:length(B)
    boundary=B{k};
    plot((boundary(:,2)-1)/100-.100,boundary(:,1)/2+5.5,'k')
end
colormap('jet')
xlim([-0.1 4.8]);ylim([6 16])
pbaspect([2,0.25,0.25]);set(gca,'TickDir','out');ax = gca;ax.YTick = frequencies(1:4:end);ax.XTick = time(11:20:end);
set(fhandle, 'Position', [100, 100, 1200, 200]);
xlabel('Time (s) - relative to memory items onset')
ylabel('Frequency (HZ)')

