%PS/NE329 Group 2 Analysis of Performance
%Edited 4/25/23

clear; %clears workspace
clc; %clears terminal

%pathing to toolboxes
addpath('/Users/sethschallies/Documents/fieldtrip-20230223');
ft_defaults


input_dataFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Preprocessed Data'; % needs updating if another computer
figFolder = '/Users/sethschallies/Desktop/PS:NE 329/Experiments/Processed Data/Figures'; % needs updating if another computer

filenames = dir(fullfile(input_dataFolder, '*.mat'));
names = string({filenames.name});
performance = zeros(209*length(filenames), 3);
i_storage = 0;
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
    for i = 1:length(output.Data)
        trial = output.Data(i);
        if ismember(i, output.EEG_stim_locked.rejected_trials) == false
            if trial.Resp_Type == 1 % if response is optimal and it is a real trial, accuracy = 1
                accuracy = 1;
                trial_num = i;
            else
                accuracy = 0;
                trial_num = i;
            end
        else % if it is not a valid trial, or if the response is none, accuracy = 3
            accuracy = NaN;
            trial_num = i;
        end
        performance(i + i_storage, 1) = sbj; % subject number
        performance(i + i_storage, 2) = i; % trial number
        performance(i + i_storage, 3) = accuracy; % binary accuracy rating
    end
    i_storage = i_storage + i;
end
%plotting

