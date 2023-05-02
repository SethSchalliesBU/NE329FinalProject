%PS329_Group2_Analysis
%last edited 4/25/23

clear;
clc;
%addpath to necessary toolboxes
addpath('/Users/sethschallies/Documents/fieldtrip-20230223');
ft_defaults

%load channel location struct
load('Chan_locs_acticap_vWMExp2a.mat')

%select plotting options
plot_indv = 1;
plot_grandAvg = 1;

input_dataFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Preprocessed Data';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE THIS PATH to where your preprocessed data is saved
figFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Processed Data/Figures';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE THIS PATH to where you want to save your figures

EOG_labels = {'TVEOG','BVEOG','LHEOG','RHEOG','StimTrak'};

%Specify Channel of Interest (choi) --> FRN is usually calculated at FCz electrode.
choi = 'FCz';

filenames = dir(fullfile(input_dataFolder, '*.mat'));
names = string({filenames.name});

%initialize variable for saving FCz ERP and FRN from each subject
StimOn_ERP = zeros(length(filenames),60,3001);%3001 is the length of the epochs in samples (or ms)
FdbOn_ERP_pos = zeros(length(filenames),60,3001);
FdbOn_ERP_neg = zeros(length(filenames),60,3001);
FRN = zeros(1,length(filenames));

%FRN Type 2. Here we consider only optimal response trials. Furthermore, we consider neutral feedback on reward trials as being
%negative feeback. On Punish trials, neutral feedback is considered as positive feedback.  
FdbOn_ERP_pos2 = zeros(length(filenames),60,3001);
FdbOn_ERP_neg2 = zeros(length(filenames),60,3001);
FRN2 = zeros(1,length(filenames));


%initialize variables for saving TF Power
TF_Pow_StimLocked = cell(length(filenames),1);
TF_Pow_FdbLocked = cell(length(filenames),1);


for sbj = 1:length(filenames)
    disp(['Sbj #' num2str(sbj)])
    load(fullfile(input_dataFolder, filenames(sbj).name));

    valid_trials = ones(1,length(output.EEG_stim_locked.sampleinfo));
    valid_trials(output.EEG_stim_locked.rejected_trials) = 0;
    Trl_Type = [output.Data.Trl_Type];%1=reward, 2=punish
    Trl_Type = Trl_Type(logical(valid_trials));
    Resp_Type = [output.Data.Resp_Type]; %1=optimal response, 2=suboptimal, 3=none
    Resp_Type = Resp_Type(logical(valid_trials));
    Fdb_Type = [output.Data.Fdb_Type]; %1=Positive, 2=negative, 3=neutral
    Fdb_Type = Fdb_Type(logical(valid_trials));
    time = output.EEG_stim_locked.time{1};

    chan_labels = output.EEG_stim_locked.label; %electrode labels
    [~,EOG_indices,~] = intersect(chan_labels,EOG_labels);
    EEG_chan_indices = 1:length(chan_labels);
    chan_labels(EOG_indices) = [];
    choi_idx = find(strcmp(chan_labels,choi)==1);
    EEG_chan_indices(EOG_indices) = [];

    %find trial indices
    pos_feedback_trial_idx = find(Fdb_Type == 1);
    neg_feedback_trial_idx = find(Fdb_Type == 2);
    neutral_feedback_trial_idx = find(Fdb_Type == 3);
   
    %find trial indices for FRN Type 2
    Reward_pos_feedback_trial_idx2 = sort([find(Trl_Type == 1 & Fdb_Type == 1 & Resp_Type == 1)]); %here we find the indices of reward trials with optimal responses and positive feedback 
    Punish_pos_feedback_trial_idx2 = sort([find(Trl_Type == 2 & Fdb_Type == 3 & Resp_Type == 1)]); %here we find the indices of punish trials with optimal responses and neutral(positive) feedback 

    Reward_neg_feedback_trial_idx2 = sort([find(Trl_Type == 1 & Fdb_Type == 3 & Resp_Type == 1)]); %here we find the indices of reward trials with optimal responses and neutral(negative) feedback
    Punish_neg_feedback_trial_idx2 = sort([find(Trl_Type == 2 & Fdb_Type == 2 & Resp_Type == 1)]); %here we find the indices of punish trials with optimal responses and negative feedback 


    baseline_timeWindow = [-.3 -.1];
    baseline_idx = dsearchn(time',baseline_timeWindow');

    %concatenate data into one matrix with dimensions num_channels x time x n_trials
    eeg_data = cat(3,output.EEG_fdb_locked.trial{:});

    FdbOn_ERP_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,pos_feedback_trial_idx),3)); %FCz erp
    FdbOn_ERP_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,neg_feedback_trial_idx),3)); %FCz erp


    FdbOn_ERP_Reward_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Reward_pos_feedback_trial_idx2),3)); 
    FdbOn_ERP_Reward_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Reward_neg_feedback_trial_idx2),3)); 
    FdbOn_ERP_Punish_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Punish_pos_feedback_trial_idx2),3)); 
    FdbOn_ERP_Punish_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Punish_neg_feedback_trial_idx2),3)); 



    %The FRN can be, and has been, computed many different ways. Here we
    %use the Difference Wave approach (Osinksy et al.,2012). We subtract
    %the ERP for postive feedback trials from the ERP from the negative
    %feedback trials. Then we take the mean amplitude of the resulting
    %difference wave between 228-344ms post feedback onset (this interval
    %was previously identified in a metaanalysis on FRN literature
    %(Sambrook and Goslin, 2015) (Mushtaq et al.,2016).
    toi = [.228 .344];
    toi_idx = dsearchn(time',toi');
    DiffWave(sbj,:,:) = FdbOn_ERP_neg(sbj,:,:) - FdbOn_ERP_pos(sbj,:,:);
    FRN(sbj) = mean(DiffWave(sbj,choi_idx,toi_idx(1):toi_idx(2)));


    DiffWave_Reward(sbj,:,:) = FdbOn_ERP_Reward_neg(sbj,:,:) - FdbOn_ERP_Reward_pos(sbj,:,:);
    FRN_Reward(sbj) = mean(DiffWave_Reward(sbj,choi_idx,toi_idx(1):toi_idx(2)));

    DiffWave_Punish(sbj,:,:) = FdbOn_ERP_Punish_neg(sbj,:,:) - FdbOn_ERP_Punish_pos(sbj,:,:);
    FRN_Punish(sbj) = mean(DiffWave_Punish(sbj,choi_idx,toi_idx(1):toi_idx(2)));


    %Accuracy
    Overall_accuracy(sbj) = sum(Resp_Type==1)/sum(Fdb_Type~=4);

    %Make sure channel order of eeg data matches channel order needed for plotting topoplots
    Chan_order_needed = {Chan_locs_acticap_vWMExp2a.labels}';
    chan_order = zeros(1,length(Chan_order_needed));
    for ch = 1:length(Chan_order_needed)
        chan_order(ch) = find(strcmp(Chan_order_needed(ch),chan_labels));
    end



    %% Time-Frequency Analysis to inspect Fontral Midline Theta (FMT) (sepertely for Stim-Locked and Fdb-Locked data)
    cfg = [];
    cfg.trials = 'all';
    cfg.channel = chan_labels;
    [eeg_data_stimLocked] = ft_selectdata(cfg, output.EEG_stim_locked);
    [eeg_data_fdbLocked] = ft_selectdata(cfg, output.EEG_fdb_locked);

    eeg_data_stimLocked = cat(3,eeg_data_stimLocked.trial{:});
    eeg_data_fdbLocked = cat(3,eeg_data_fdbLocked.trial{:});

    [nchan,npnts,ntrials] = size(eeg_data_stimLocked);
    freq_range = [3, 54];

    %wavelet params
    srate = output.EEG_stim_locked.fsample;%sampling frequency
    wavtime = -2:1/srate:2;
    t = npnts;
    nData = (t)*ntrials;
    nKern = length(wavtime);
    nConv = nData+nKern-1;
    halfwav = floor(length(wavtime)/2)+1;
    nFrex = 30;
    frex = exp(linspace(log(freq_range(1)),log(freq_range(2)),nFrex));
    s = linspace(4,12,nFrex)./ (2*pi.*frex);

    fft_data_StimLocked = fft(reshape(eeg_data_stimLocked,nchan,t*ntrials),nConv,2);
    fft_data_FdbLocked = fft(reshape(eeg_data_fdbLocked,nchan,t*ntrials),nConv,2);

    bsl_win_TF = [-.4, -.2]; %baseline time window for TF analysis

    bsl_st_idx = dsearchn(time',bsl_win_TF(1));%idx of baseline window start
    bsl_end_idx = dsearchn(time',bsl_win_TF(1));%idx of baseline window end

    TF_Power_StimLocked = zeros(nFrex,npnts);
    TF_Power_FdbLocked = zeros(nFrex,npnts);


    for fi = 1:nFrex
        cmw = exp(2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2));
        cmwX = fft(cmw,nConv);
        cmwX = cmwX./max(cmwX);

        %compute baseline normalized Total Power for StimLocked Data at FCz
        as = ifft(fft_data_StimLocked(choi_idx,:) .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_StimLocked = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_StimLocked(bsl_st_idx:bsl_end_idx));
        %         TF_Power_StimLocked(fi,:) = 100*((tmp_TF_Power_StimLocked-bsl)/bsl);% percent change from baseline
        TF_Power_StimLocked(fi,:) = 10*log10(tmp_TF_Power_StimLocked./bsl);% dB conversion

        as = ifft(fft_data_FdbLocked(choi_idx,:) .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_FdbLocked = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_FdbLocked(bsl_st_idx:bsl_end_idx));
        %         TF_Power_FdbLocked(fi,:) = 100*((tmp_TF_Power_FdbLocked-bsl)/bsl);% percent change from baseline
        TF_Power_FdbLocked(fi,:) = 10*log10(tmp_TF_Power_FdbLocked./bsl);% dB conversion
    end
    TF_Pow_StimLocked{sbj} = TF_Power_StimLocked;
    TF_Pow_FdbLocked{sbj} = TF_Power_FdbLocked;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plotting individual subjects

    if plot_indv
        %plotting difference wave
        fig=figure('Position', [80 80 800 400]);
        tl = tiledlayout(1,2);
        tl.TileSpacing = 'compact';
        title(tl,['Sbj # ' num2str(sbj)])

        nexttile
        plot(time,squeeze(FdbOn_ERP_pos(sbj,choi_idx,:)),'b');hold on;
        plot(time,squeeze(FdbOn_ERP_neg(sbj,choi_idx,:)),'r');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title(' Pos Fdb ERP vs Neg Fdb ERP')

        nexttile
        plot(time,squeeze(DiffWave(sbj,choi_idx,:)),'b');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title(' Difference Wave')

        mkdir(fullfile(figFolder,'Individual Sbj Plots'))
        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave.jpg']));
        close(fig)

        %plotting difference wave For Reward Trials 
        fig=figure('Position', [80 80 800 400]);
        tl = tiledlayout(1,2);
        tl.TileSpacing = 'compact';
        title(tl,['Sbj # ' num2str(sbj)])

        nexttile
        plot(time,squeeze(FdbOn_ERP_Reward_pos(sbj,choi_idx,:)),'b');hold on;
        plot(time,squeeze(FdbOn_ERP_Reward_neg(sbj,choi_idx,:)),'r');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({' Pos Fdb ERP vs Neg Fdb ERP', 'Reward Trials'})

        nexttile
        plot(time,squeeze(DiffWave_Reward(sbj,choi_idx,:)),'b');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({'Difference Wave', 'Reward Trials'})

        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave_Reward.jpg']));
        close(fig)

        %plotting difference wave For Punish Trials 
        fig=figure('Position', [80 80 800 400]);
        tl = tiledlayout(1,2);
        tl.TileSpacing = 'compact';
        title(tl,['Sbj # ' num2str(sbj)])

        nexttile
        plot(time,squeeze(FdbOn_ERP_Punish_pos(sbj,choi_idx,:)),'b');hold on;
        plot(time,squeeze(FdbOn_ERP_Punish_neg(sbj,choi_idx,:)),'r');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({' Pos Fdb ERP vs Neg Fdb ERP', 'Punish Trials'})

        nexttile
        plot(time,squeeze(DiffWave_Punish(sbj,choi_idx,:)),'b');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({'Difference Wave', 'Punish Trials'})

        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave_Punish.jpg']));
        close(fig)


        %plot topoplots for Pos fdb Neg fdb and DiffWave
        fig=figure('Position', [80 80 800 400],'Name',['Sbj' num2str(sbj)]);
        tl = tiledlayout(1,3);
        tl.TileSpacing = 'compact';

        nexttile
        temp_data = squeeze(mean(FdbOn_ERP_neg(sbj,:,toi_idx(1):toi_idx(2)),3));
        topoplot(temp_data',Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
        colormap('turbo')
        c = colorbar;
        c.Label.String = '\mu V';
        title(['Neg Fdb'])
        caxis([-5 5]);

        nexttile
        temp_data = squeeze(mean(FdbOn_ERP_pos(sbj,:,toi_idx(1):toi_idx(2)),3));
        topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
        colormap('turbo')
        c = colorbar;
        c.Label.String = '\mu V';
        title(['Pos Fdb'])
        caxis([-5 5]);

        nexttile
        temp_data = squeeze(mean(DiffWave(sbj,:,toi_idx(1):toi_idx(2)),3));
        topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
        colormap('turbo')
        c = colorbar;
        c.Label.String = '\mu V';
        title(['Neg-Pos'])
        caxis([-5 5]);

        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_Topoplots.jpg']));
        close(fig)

        %TF-plot for visualizeing FMT
        fig=figure;
        contourf(time,frex,TF_Power_FdbLocked,40,'linecolor','none');
        set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
        colormap(turbo)
        ylabel('Freq [Hz]')
        xlabel('Time from Fdb Onset')
        xlim([-.5 2])
        caxis([-5 5])
        a=colorbar;
        ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
        a.Label.Position(1) = 3;
        a.Ticks = [-5,0,5];
        title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})

        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_FdbLocked.jpg']));
        close(fig)


        fig=figure;
        contourf(time,frex,TF_Power_StimLocked,40,'linecolor','none');
        set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
        colormap(turbo)
        ylabel('Freq [Hz]')
        xlabel('Time from Stim Onset')
        xlim([-.5 2])
        caxis([-5 5])
        a=colorbar;
        ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
        a.Label.Position(1) = 3;
        a.Ticks = [-5,0,5];
        title({['Sbj #' num2str(sbj) ' TF-Power'], 'Stimulus Locked - FCz'})

        saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_StimLocked.jpg']));
        close(fig)
    end
end

if plot_grandAvg
    %% Plotting Group Average

    %plotting difference wave
    fig=figure('Position', [80 80 800 400]);
    tl = tiledlayout(1,2);
    tl.TileSpacing = 'compact';
    title(tl,'Grand Average')

    nexttile
    plot(time,squeeze(mean(FdbOn_ERP_pos(:,choi_idx,:),1)),'b');hold on;
    plot(time,squeeze(mean(FdbOn_ERP_neg(:,choi_idx,:),1)),'r');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Pos Fdb ERP vs Neg Fdb ERP')

    nexttile
    plot(time,squeeze(mean(DiffWave(:,choi_idx,:),1)),'b');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Difference Wave')

    mkdir(fullfile(figFolder,'Grand Average'))
    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave.jpg'));
    close(fig)

    %plotting difference Reward Trials
    fig=figure('Position', [80 80 800 400]);
    tl = tiledlayout(1,2);
    tl.TileSpacing = 'compact';
    title(tl,'Grand Average')

    nexttile
    plot(time,squeeze(mean(FdbOn_ERP_Reward_pos(:,choi_idx,:),1)),'b');hold on;
    plot(time,squeeze(mean(FdbOn_ERP_Reward_neg(:,choi_idx,:),1)),'r');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title({'Pos Fdb ERP vs Neg Fdb ERP', 'Reward Trials'})

    nexttile
    plot(time,squeeze(mean(DiffWave_Reward(:,choi_idx,:),1)),'b');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title({'Difference Wave', 'Reward Trials'})

    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave_Reward.jpg'));
    close(fig)

    %plotting difference Punish Trials
    fig=figure('Position', [80 80 800 400]);
    tl = tiledlayout(1,2);
    tl.TileSpacing = 'compact';
    title(tl,'Grand Average')

    nexttile
    plot(time,squeeze(mean(FdbOn_ERP_Punish_pos(:,choi_idx,:),1)),'b');hold on;
    plot(time,squeeze(mean(FdbOn_ERP_Punish_neg(:,choi_idx,:),1)),'r');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title({'Pos Fdb ERP vs Neg Fdb ERP', 'Punish Trials'})

    nexttile
    plot(time,squeeze(mean(DiffWave_Punish(:,choi_idx,:),1)),'b');
    xlim([-.2 .5])
    ylabel('Amplitude [\mu V]')
    xlabel('Time post Fdb [ms]')
    x = [toi(1) toi(1) toi(2) toi(2)];
    y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
    b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    X_zeroLine = xline(0,'k');
    X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Y_zeroLine = yline(0,'k');
    Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title({'Difference Wave', 'Punish Trials'})

    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave_Punish.jpg'));
    close(fig)


    %plot topoplots for Pos fdb Neg fdb and DiffWave
    fig=figure('Position', [80 80 800 400],'Name','Grand Average');
    tl = tiledlayout(1,3);
    tl.TileSpacing = 'compact';

    nexttile
    temp_data = squeeze(mean(mean(FdbOn_ERP_neg(:,:,toi_idx(1):toi_idx(2)),3),1));
    topoplot(temp_data',Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
    colormap('turbo')
    c = colorbar;
    c.Label.String = '\mu V';
    title(['Neg Fdb'])
    caxis([-5 5]);

    nexttile
    temp_data = squeeze(mean(mean(FdbOn_ERP_pos(:,:,toi_idx(1):toi_idx(2)),3),1));
    topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
    colormap('turbo')
    c = colorbar;
    c.Label.String = '\mu V';
    title(['Pos Fdb'])
    caxis([-5 5]);

    nexttile
    temp_data = squeeze(mean(mean(DiffWave(:,:,toi_idx(1):toi_idx(2)),3),1));
    topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
    colormap('turbo')
    c = colorbar;
    c.Label.String = '\mu V';
    title(['Neg-Pos'])
    caxis([-5 5]);

    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_Topoplots.jpg'));
    close(fig)

    %TF-plot for visualizeing FMT
    fig=figure;
    contourf(time,frex,squeeze(mean(cat(3,TF_Pow_FdbLocked{:}),3)),40,'linecolor','none');
    set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
    colormap(turbo)
    ylabel('Freq [Hz]')
    xlabel('Time from Fdb Onset')
    xlim([-.5 2])
    caxis([-5 5])
    a=colorbar;
    ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
    a.Label.Position(1) = 3;
    a.Ticks = [-5,0,5];
    title({'Grand Average TF-Power', 'Feedback Locked - FCz'})

    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_FdbLocked.jpg'));
    close(fig)


    fig=figure;
    contourf(time,frex,squeeze(mean(cat(3,TF_Pow_StimLocked{:}),3)),40,'linecolor','none');
    set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','lo g','YMinorTick','off');
    colormap(turbo)
    ylabel('Freq [Hz]')
    xlabel('Time from Stim Onset')
    xlim([-.5 2])
    caxis([-5 5])
    a=colorbar;
    ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
    a.Label.Position(1) = 3;
    a.Ticks = [-5,0,5];
    title({'Grand Average TF-Power', 'Stimulus Locked - FCz'})

    saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_StimLocked.jpg'));
    close(fig)
end

%% Stats?
