function [clean_eeg_data, Trial_Info] = EEG_PreProcessing_Group2(eeg_filename,dataFolder,automatic_trl_rej)
%%EEG_PreProcessing_Group2
% By CTGill
% Last updated 4/18/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Trials
tic
cfg = [];
cfg.dataset = [dataFolder filesep eeg_filename];
cfg.trialfun = 'PRLTask_trialfun';
evalc('TrlDefs = ft_definetrial(cfg)'); 
Trial_Info = TrlDefs.trl(:,[4:10]);
fullcfg = rmfield(TrlDefs,'trl');
evalc('data_eeg = ft_preprocessing(fullcfg)');
display(['Defining trials took ' num2str(toc) ' seconds to complete'])

%% detrending and filtering
tic
cfg = [];
cfg.demean = 'yes';
cfg.bpfilter = 'yes';       % enable band-pass filtering
cfg.bpfilttype = 'firws';
cfg.bpfreq = [1 60];        %[High-pass freq, Low-pass freq];
cfg.detrend = 'yes';
evalc('data_eeg = ft_preprocessing(cfg,data_eeg)');

display(['Detrending and filering dataset took ' num2str(toc) ' seconds to complete'])

%% Apply trial definition to data.
tic
Trialcfg.trl = TrlDefs.trl;
evalc('data_eeg = ft_redefinetrial(Trialcfg,data_eeg)');

display(['Applying trial definition to data took ' num2str(toc) ' seconds to complete'])

%% Trial Rejection\Channel Rejection
if automatic_trl_rej == 0

    %Manual
    %manually inspect every trial marking any artifacts. If there are any
    %noisy channels that, by themselves, would cause 5-10% of all trials to be
    %removed, then keep these trials and remove that channel. 
    cfg = [];
    cfg.ylim = [0 20];
    cfg.preproc.bpfilter = 'yes';
    cfg.preproc.bpfreq = [4 60];
    cfg.channel = data_eeg.label;  
    cfg.viewmode = 'vertical';
    cfg.allowoverlap = 'yes';
    data_rej = ft_databrowser(cfg,data_eeg);

    %find any trials with artifacts and remove these from further analysis
    trials2remove = [];
    for art = 1:size(data_rej.artfctdef.visual.artifact,1)
        trials2remove(art) = find(data_eeg.sampleinfo(:,1)<= data_rej.artfctdef.visual.artifact(art,1) & data_eeg.sampleinfo(:,2)>= data_rej.artfctdef.visual.artifact(art,1));
    end

    channels2remove = [];
    prompt = {'Enter index of any channels that should be removed separated by a space'}; %description of fields
    prompt_title = 'Bad Channels to Remove';
    dims = [1 50];
    defaults = {' '};
    answer = inputdlg(prompt,prompt_title,dims,defaults);
    channels2remove = str2num(answer{1,:}); %Gets additional ICs to remove


else
    %%Automatic
    %excluded trials if the peak-to-peak voltage within the EEG epoch was greater than
    %300uV in any 200-ms window in any channel. Bacigalupo & Luck (2020) excluded
    %trials if the peak-to-peak voltage within the EEG epoch was greater than
    %300uV in any 200-ms window in any channel. 
    
    tic
    art=ones(1,400);
    for trl = 1:400
        for i = 1:size(data_eeg.trial{trl},2)-200
            p2p = range(data_eeg.trial{trl}(:,i:i+200),2);
            if numel(find(p2p>=300))>0
                art(trl)=0;
                break;
            end
        end
    end
    channels2remove = [];
    trials2remove = find(art==0);
end

data_eeg.trialinfo(trials2remove,:)=[];
data_eeg.sampleinfo(trials2remove,:)=[];
data_eeg.trial(:,trials2remove)=[];
data_eeg.time(:,trials2remove)=[];
data_eeg.rejected_trials = trials2remove;
data_eeg.rejected_channels = channels2remove;

display(' ')
display('---Trial Rejection Applied---')
display(['Rejected ' num2str(length(trials2remove)) ' trials'])
display(['Rejected Trial #s ' num2str([trials2remove(:)'])])
display(['Remaining ' num2str(length(data_eeg.trial)) ' trials'])
display(['Remaining ' num2str((length(data_eeg.trial)/(size(data_eeg.trialinfo,1)+numel(data_eeg.rejected_trials)))*100) '% trials'])
display(' ')
display(['Trial Rejection took ' num2str(toc) 'seconds to complete'])

%% Run ICA
tic
cfg_ica = [];
cfg_ica.method = 'runica';
data_eeg_ica = rmfield(ft_componentanalysis(cfg_ica, data_eeg),'cfg');

display(['ICA took ' num2str(toc) ' seconds to complete'])
%% Selecting ICs to remove
%manually finding additional ICA Components that should be removed. The most commonly 
% removed components are blink components and components contain vertical and horizontal eye movements.
% If present, heart beat component and bad channel components should also be removed.
% Only need to inspect the first 20 componenets because they account for a
% majority of the signal variance (see mike x cohen
%https://mikexcohen.com/lecturelets/ica2clean/ica2clean.html)
%
%1)look at topoplots and the time course 1st to find potential components
%to remove and mark these for further inspection. 
%
%2) look at ERP image to see if there is task related signal that is
%present. Also look at activity power spectrogram to see if it looks like
%the normal 1/f distribution.
%
%3) after inspecting the ERP image and power spectrum for each IC that you
%are considering removing, input 'y' into the command window if the component
%being inspected should be removed or 'n' if the component should not be
%removed and then press 'enter' to show the plots for the next IC marked
%for inspection.

cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'component';
if strcmp(data_eeg.label{1},'FTT9h')
    cfg.elec = ft_read_sens('standard_1020_CG_modified.elc');
else
    cfg.elec = ft_read_sens('standard_1020.elc');
end
cfg.allowoverlap = 'yes';
component_rej = ft_databrowser(cfg,data_eeg_ica);

prompt = {'Enter ICs that may need to be removed separated by a space'};
title1 = 'ICs to Remove';
dims = [1 50];
defaults = {' '};
answer = inputdlg(prompt,title1,dims,defaults);
IC2remove = str2num(answer{1,:}); %Gets ICs that need to be remove and ICs that require further inspection 

clearvars ERPimages
ERPimages = cell(1,numel(IC2remove));
%%Generate ERPimage for each component in question after locking to feedback onset
for IC = 1:numel(IC2remove)
    %find the min epoch length in case they are not all exactly the same
    min_epoch_length = [];
    for trl = 1:length(data_eeg_ica.trial)
        if isempty(min_epoch_length) || (length(data_eeg_ica.trial{trl}) < min_epoch_length)
            min_epoch_length = length(data_eeg_ica.trial{trl});
        end
    end

    tmp_comp=[];
    for i = 1:length((data_eeg_ica.trial))
        tmp_comp=[tmp_comp;data_eeg_ica.trial{i}(IC2remove(IC),1:min_epoch_length)];
    end
    ERPimages{IC} = tmp_comp;
end

%%generate power spectrum for each IC being inspected
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = IC2remove;%'all';%compute the power spectrum in all ICs
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 0:1:60;
freq = ft_freqanalysis(cfg, data_eeg_ica);

%%plot
for IC = 1: length(ERPimages)
    %plot ERPimage for each IC in question
    fig = figure('Position', [50 50 500 750],'Name',['IC #' num2str(IC2remove(IC))]);
    t = tiledlayout(6,1);
    t.TileSpacing = 'compact';
    t.Title.String = ['IC #' num2str(IC2remove(IC))];

    nexttile;
    plot([data_eeg_ica.trial{1}(IC2remove(IC),:),data_eeg_ica.trial{2}(IC2remove(IC),:),data_eeg_ica.trial{3}(IC2remove(IC),:)])
    xlabel('time [ms]');
    ylabel('Amp [\muV]')
    set(gca,'xticklabel',{[]});
    title('IC timecourse: 1st three trials')

    nexttile(2,[2,1]);
    imagesc(ERPimages{IC});
    xticks(500:500:3500)
    xticklabels({'','0','','','','',''})
    c=colorbar;
    xlabel('trial #');
    ylabel('time [ms]');

    nexttile;
    plot(mean(ERPimages{IC}))
    xline(1000,'--r')
    xticks(500:500:3500)
    xticklabels({'','0','','','','',''})
    title('IC ERP')
    
    nexttile([2,1]);
    cfg.channel = IC;
    ft_singleplotER(cfg,freq);
    plot(freq.freq,10*log10(freq.powspctrm(IC,:)))
    xlabel('Freq [Hz]');
    ylabel('Power 10*log10');
    xlim([1 58]);


    figure;
    cfg = [];
    cfg.component = IC2remove(IC);       % specify the component(s) that should be
    if strcmp(data_eeg.label{1},'FTT9h')
        cfg.elec = ft_read_sens('standard_1020_CG_modified.elc');
    else
        cfg.elec = ft_read_sens('standard_1020.elc');
    end
    cfg.comment   = 'no';
    cfg.colorbar  = 'yes';
    ft_topoplotIC(cfg, data_eeg_ica);
   
    pause(1)
    %continue to next trial or break
    a = input('Should this IC be removed (y/n)? ','s');
    if strcmpi(a,'y')
        close all;
    else
        IC2remove(IC) = nan;
        close all;
    end
end

IC2remove = IC2remove(~isnan(IC2remove));

%%remove the bad components and backproject the data
cfg = [];
cfg.component = IC2remove; % to be removed component(s)
clean_eeg_data = ft_rejectcomponent(cfg, data_eeg_ica);

%% Interpolate removed electrodes
if ~isempty(channels2remove)
    disp(['interpolating over ' num2str(numel(channels2remove)) ' removed channels'])

    cfg_neighb          = [];
    cfg_neighb.feedback = 'no';
    cfg_neighb.method   = 'triangulation';
    if strcmp(data_eeg.label{1},'FTT9h')
    cfg_neighb.elec     = ft_read_sens('standard_1020_CG_modified.elc');
    else
    cfg_neighb.elec     = ft_read_sens('standard_1020.elc');
    end
    neighbours          = ft_prepare_neighbours(cfg_neighb);

    cfg=[];
    cfg.channel = [clean_eeg_data.label];
    temp_data = ft_selectdata(cfg,clean_eeg_data);
    
    cfg = [];
    cfg.method = 'spline';
    if strcmp(data_eeg.label{1},'FTT9h')
        cfg.elec = ft_read_sens('standard_1020_CG_modified.elc');
    else
        cfg.elec = ft_read_sens('standard_1020.elc');
    end
    cfg.missingchannel = channels2remove;
    cfg.neighbours = neighbours;
    cfg.senstype = 'eeg';
    
    clean_eeg_data = ft_channelrepair_oldVersion(cfg,temp_data);
end
%% Implicit Rereferencing
tic
Imp_cfg = [];
Imp_cfg.reref = 'yes';
Imp_cfg.refmethod = 'avg';
Imp_cfg.channel  = 'all';
Imp_cfg.implicitref = 'TP10';
Imp_cfg.refchannel = {'TP9', 'TP10'}; %Average over the mastoids % 'all';
evalc('clean_eeg_data = ft_preprocessing(Imp_cfg,clean_eeg_data)');

display(['Rereferencing took ' num2str(toc) ' seconds to complete'])

%% Include which trials,channels, and components were removed during preprocessing
clean_eeg_data.rejected_ICs = IC2remove;
clean_eeg_data.rejected_trials = trials2remove;
clean_eeg_data.rejected_channels = channels2remove;









