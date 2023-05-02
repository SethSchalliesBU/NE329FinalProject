%% PrepareTrials_Group2
% By CTGill
% Last updated 4/18/23
%
% Before running this script you must update the paths for the
% Raw Data Folder and the Preprocessed Data Folder
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up
clear;
clc;

%add path to fieldtrip
addpath('/Users/sethschallies/Documents/fieldtrip-20230223');
ft_defaults

%Folder from which to read in Raw data
dataFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Raw Data';

%Folder for saving Preprocessed eeg data
Preprocessed_dataFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Preprocessed Data';

%set parameters
skip_already_processed_datafiles = 1; %set this value equal to 1 if you want to skip over subjects whos data has already been preprocessed and saved
automatic_trl_rej = 0;% 0=Manual trial rejection (recommended); 1=automatic trial rejection (if using this, you should redo the Trial Rejection section of EEG_PreProcessing.m to fit your desired trial exclusion criteria)


% read in filenames of unprocessed data
filenames.eeg = dir(fullfile(dataFolder, '*.eeg'));


%% Create data trial structure for each subject
for sbj=1:length(filenames.eeg)

    if skip_already_processed_datafiles
        if exist(fullfile(Preprocessed_dataFolder, [filenames.eeg(sbj).name(1:end-4) '.mat']),'file') == 2
            continue;
        end
    end

    [EEG,Trial_Info] = EEG_PreProcessing_Group2(filenames.eeg(sbj).name,dataFolder,automatic_trl_rej);


    offset = EEG.trialinfo(:,7);
    prestim = 1;
    poststim1 = 3;
    poststim2 = 2;

    cfg=[];
    cfg.begsample = 1;
    cfg.endsample = round(1 + poststim1*EEG.fsample);
    stimOn_locked_data = ft_redefinetrial(cfg,EEG);

    cfg=[];
    cfg.begsample = round(offset - prestim*EEG.fsample);
    cfg.endsample = round(offset + poststim2*EEG.fsample);
    fdbOn_locked_data = ft_redefinetrial(cfg,EEG);

    output.EEG_stim_locked = stimOn_locked_data;
    output.EEG_fdb_locked = fdbOn_locked_data;


    % Make a Trial Info Struct
    %     column 1 = trial number
    %     column 2 = trial type  (reward=1, punishment=3)
    %     column 3 = response type (1=Optimal response , 2=Suboptimal response, 3=No response)
    %     column 4 = feedback type (1=Postive Feedback, 2 = Negative Feedback, 3=Neutral Feedback, 4=if the subject did not respond in this trial)
    Data = struct('Trl_num',[],'Trl_Type',[],'Resp_Type',[],'Fdb_Type',[]);
    for trl = 1:length(Trial_Info)
        Data(trl).Trl_num = Trial_Info(trl,1);
        Data(trl).Trl_Type = Trial_Info(trl,2);
        Data(trl).Resp_Type = Trial_Info(trl,3);
        Data(trl).Fdb_Type = Trial_Info(trl,4);
    end
    output.Data = Data;


    % save output for each subject
    save([Preprocessed_dataFolder filesep [filenames.eeg(sbj).name(1:end-4) '.mat']],'output','-v7.3','-nocompression');
end