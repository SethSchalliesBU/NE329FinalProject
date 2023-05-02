%PS329_Group2_Analysis
%last edited 5/2/23
%credit: Christopher Gill, Alicia Qin, reviewed by Seth Schallies


clear;
clc;
%addpath to necessary toolboxes
addpath('C:\Users\Dell-XPS9360\Documents\PS 329\fieldtrip-20230223');
ft_defaults
addpath('C:\Users\Dell-XPS9360\Documents\PS 329\fieldtrip-20230223\eeglab_current\eeglab2022.1');

%load channel location struct
load('Chan_locs_acticap_vWMExp2a.mat')

%select plotting options
plot_indv = 1;
plot_grandAvg = 1;

input_dataFolder = 'C:\Users\Dell-XPS9360\Documents\PS 329\PS329_Experiments\Results';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE THIS PATH to where your preprocessed data is saved
figFolder = 'C:\Users\Dell-XPS9360\Documents\PS 329\PS329_Experiments\Figures1';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE THIS PATH to where you want to save your figures

EOG_labels = {'TVEOG','BVEOG','LHEOG','RHEOG','StimTrak'};

%Specify Channel of Interest (choi) --> FRN is usually calculated at FCz electrode.
choi = 'FCz';

filenames = dir(fullfile(input_dataFolder, '*.mat'));
names = string({filenames.name});

%subject and size variables
num_sbj = length(filenames);

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

TF_Pow_Reward_unexp = cell(length(filenames),1);       %add for expectation comparisons
TF_Pow_Reward_exp = cell(length(filenames),1);
TF_Pow_Punish_unexp = cell(length(filenames),1);
TF_Pow_Punish_exp = cell(length(filenames),1);

TF_Pow_Begin = cell(length(filenames),1);
TF_Pow_Mid = cell(length(filenames),1);
TF_Pow_Fin = cell(length(filenames),1);

TF_Pow_Block = {};
TF_Pow_Blk1 = cell(length(filenames),1);
TF_Pow_Blk2 = cell(length(filenames),1);
TF_Pow_Flk3 = cell(length(filenames),1);


%init variables for calculated FRN
FRN_Punish = NaN(1,num_sbj);
FRN_Reward = NaN(1,num_sbj);

FRN_Blockbegin = NaN(1,num_sbj);
FRN_Blockmid =  NaN(1,num_sbj);
FRN_Blockfin = NaN(1,num_sbj);

FRN_allbegin = NaN(1,num_sbj);
FRN_allmid = NaN(1,num_sbj);
FRN_allfin = NaN(1,num_sbj);

    %over time
    num_sec = 16;
    FRN_blocksec = cell(num_sbj,num_sec);

%init ERP (all channel) variable
FdbOn_ERP_allbegin_unexp = NaN(num_sbj,60,3001);

%initialize variable for running accuracy
%max_sec = ceil(length(all_idx)/2);
%accuracy_running = NaN(8,320);     %max 320 trials


for sbj = 1:length(filenames)
    disp(['Sbj #' num2str(sbj)])
    load(fullfile(input_dataFolder, filenames(sbj).name));

    %define valid trials and trial types
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

    tic 
    %find trial indices
    pos_feedback_trial_idx = find(Fdb_Type == 1);
    neg_feedback_trial_idx = find(Fdb_Type == 2);
    neutral_feedback_trial_idx = find(Fdb_Type == 3);
    
    %find trial indices for FRN Type 2
    Reward_expected_trial_idx = sort([find(Trl_Type == 1 & Fdb_Type == 1 & Resp_Type == 1)]); %here we find the indices of reward trials with optimal responses and positive feedback 
    Punish_expected_trial_idx = sort([find(Trl_Type == 2 & Fdb_Type == 3 & Resp_Type == 1)]); %here we find the indices of punish trials with optimal responses and neutral(positive) feedback 

    Reward_unexpected_trial_idx = sort([find(Trl_Type == 1 & Fdb_Type == 3 & Resp_Type == 1)]); %here we find the indices of reward trials with optimal responses and neutral(negative) feedback
    Punish_unexpected_trial_idx = sort([find(Trl_Type == 2 & Fdb_Type == 2 & Resp_Type == 1)]); %here we find the indices of punish trials with optimal responses and negative feedback 

    Expected_all_trial_idx = cat(2, Reward_expected_trial_idx, Punish_expected_trial_idx);
    Unexpected_all_trial_idx = cat(2, Reward_unexpected_trial_idx, Punish_unexpected_trial_idx);
    
    %all trials --> for computing across time
    all_idx = sort(cat(2, Expected_all_trial_idx, Unexpected_all_trial_idx));


    %TIME: block-3rds
    %BLOCKS
    block_len = 80;
    num_blocks = ceil(max(all_idx)/block_len);

    %create matrix of block indices (block, trials)
    block_idx = zeros(num_blocks, block_len);
    block_start = 1;
    block_end = block_start + block_len -1;
    for b = 1:num_blocks   %each block
        tr1 = 1;
        for i = 1:length(all_idx) 
                if (block_start <= all_idx(i) && all_idx(i) <= block_end && block_idx(b,tr1) == 0)
                    block_idx(b,tr1) = all_idx(i);
                    %fprintf("--added %d to (%d,%d) \n", all_idx(i),b,tr1)
                    tr1 = tr1 + 1;
                end
        end
        block_start = block_start + block_len;
        block_end = block_end + block_len;
        %fprintf("****going to next block: %d \n", block_start)
    end
    

    %divide into begin/mid/final trials for each block
    third_block = ceil(block_len/3);                         
    block_begin = zeros();
    block_mid = zeros();
    block_fin = zeros();

    %-->splitting by content
    for i= 1:num_blocks
        block_beg_start = block_len*(i-1);
        block_mid_start = block_len*(i-1) + third_block + 1;
        block_fin_start = block_len*(i-1) + 2*third_block + 1;
        block_fin_end = block_fin_start + third_block - 1;
        for j = 1:block_len
            if (block_beg_start <= block_idx(i,j) && block_idx(i,j) < block_mid_start)
                block_begin = cat(2, block_begin, block_idx(i,j));
            elseif (block_mid_start <= block_idx(i,j) && block_idx(i,j)< block_fin_start)
                block_mid = cat(2, block_mid, block_idx(i,j));
            elseif (block_fin_start <= block_idx(i,j) && block_idx(i,j) <= block_fin_end)
                block_fin = cat(2, block_fin, block_idx(i,j));
            end
        end
    end

%     %convert to vector
%     block_begin = sort(reshape(block_begin, 1, []));
%     block_mid = sort(reshape(block_mid, 1, []));
%     block_fin = sort(reshape(block_fin, 1, []));

    %remove zeros
    block_begin = find(block_begin ~= 0);
    block_mid = find(block_mid ~= 0);
    block_fin = find(block_fin ~= 0);
    

    %block section & expectation trial indices
    block_begin_punish_unexp_idx = intersect(block_begin, Punish_unexpected_trial_idx)';
    block_begin_punish_exp_idx = intersect(block_begin, Punish_expected_trial_idx)';
    block_begin_reward_unexp_idx = intersect(block_begin, Reward_unexpected_trial_idx)';
    block_begin_reward_exp_idx = intersect(block_begin, Reward_expected_trial_idx)';

    block_mid_punish_unexp_idx = intersect(block_mid, Punish_unexpected_trial_idx)';
    block_mid_punish_exp_idx = intersect(block_mid, Punish_expected_trial_idx)';
    block_mid_reward_unexp_idx = intersect(block_mid, Reward_unexpected_trial_idx)';
    block_mid_reward_exp_idx = intersect(block_mid, Reward_expected_trial_idx)';

    block_fin_punish_unexp_idx = intersect(block_fin, Punish_unexpected_trial_idx)';
    block_fin_punish_exp_idx = intersect(block_fin, Punish_expected_trial_idx)';
    block_fin_reward_unexp_idx = intersect(block_fin, Reward_unexpected_trial_idx)';
    block_fin_reward_exp_idx = intersect(block_fin, Reward_expected_trial_idx)';

    %TIME: finer grain across blocks -- 
    time_idx = {};
        %divide into begin/mid/final trials for each block
        fifth_block = ceil(block_len/5);      %16 trials per section
        block_sec1 = zeros();
        block_sec2 = zeros();
        block_sec3 = zeros();
        block_sec4 = zeros();
        block_sec5 = zeros();

        %-->splitting by content
        for i= 1:num_blocks
            block_beg_start = block_len*(i-1);
            block_sec1_start = block_len*(i-1) + fifth_block + 1;
            block_sec2_start = block_len*(i-1) + 2*fifth_block + 1;
            block_sec3_start = block_len*(i-1) + 3*fifth_block + 1;
            block_sec4_start = block_len*(i-1) + 4*fifth_block + 1;
            block_sec5_start = block_len*(i-1) + 5*fifth_block + 1;
            block_sec5_end = block_sec5_start + fifth_block - 1;
            for j = 1:block_len
                if (block_beg_start <= block_idx(i,j) && block_idx(i,j) < block_sec2_start)
                    block_sec1 = cat(2, block_sec1, block_idx(i,j));
                elseif (block_sec2_start <= block_idx(i,j) && block_idx(i,j)< block_sec3_start)
                    block_sec2 = cat(2, block_sec2, block_idx(i,j));
                elseif (block_sec3_start <= block_idx(i,j) && block_idx(i,j) <= block_sec4_start)
                    block_sec3 = cat(2, block_sec3, block_idx(i,j));
                elseif (block_sec4_start <= block_idx(i,j) && block_idx(i,j) <= block_sec5_start)
                    block_sec4 = cat(2, block_sec4, block_idx(i,j));
                elseif (block_sec5_start <= block_idx(i,j) && block_idx(i,j) <= block_sec5_end)
                    block_sec5 = cat(2, block_sec5, block_idx(i,j));
                end
            end
        end
    
        %remove zeros
        block_sec1 = find(block_sec1 ~= 0);
        block_sec2 = find(block_sec2 ~= 0);
        block_sec3 = find(block_sec3 ~= 0);
        block_sec4 = find(block_sec4 ~= 0);
        block_sec5 = find(block_sec5 ~= 0);



 %% ___
    baseline_timeWindow = [-.3 -.1];
    baseline_idx = dsearchn(time',baseline_timeWindow');


    %concatenate data into one matrix with dimensions num_channels x time x n_trials
    eeg_data = cat(3,output.EEG_fdb_locked.trial{:});

    FdbOn_ERP_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,pos_feedback_trial_idx),3)); %FCz erp
    FdbOn_ERP_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,neg_feedback_trial_idx),3)); %FCz erp


    FdbOn_ERP_Reward_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Reward_expected_trial_idx),3)); 
    FdbOn_ERP_Reward_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Reward_unexpected_trial_idx),3)); 
    FdbOn_ERP_Punish_pos(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Punish_expected_trial_idx),3)); 
    FdbOn_ERP_Punish_neg(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,Punish_unexpected_trial_idx),3)); 


    %ERP by time section within block

    FdbOn_ERP_blockbegin_punish_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_begin_punish_unexp_idx),3)); 
    FdbOn_ERP_blockbegin_punish_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_begin_punish_exp_idx),3)); 
    FdbOn_ERP_blockbegin_reward_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_begin_reward_unexp_idx),3)); 
    FdbOn_ERP_blockbegin_reward_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_begin_reward_exp_idx),3)); 

    FdbOn_ERP_blockmid_punish_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_mid_punish_unexp_idx),3)); 
    FdbOn_ERP_blockmid_punish_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_mid_punish_exp_idx),3)); 
    FdbOn_ERP_blockmid_reward_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_mid_reward_unexp_idx),3)); 
    FdbOn_ERP_blockmid_reward_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_mid_reward_exp_idx),3)); 

    FdbOn_ERP_blockfin_punish_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_fin_punish_unexp_idx),3)); 
    FdbOn_ERP_blockfin_punish_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_fin_punish_exp_idx),3)); 
    FdbOn_ERP_blockfin_reward_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_fin_reward_unexp_idx),3)); 
    FdbOn_ERP_blockfin_reward_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,block_fin_reward_exp_idx),3)); 

    % 3 blocksec X 2 trialtype X 2 fdb_expectation
    FdbOn_ERP_blocksecXtrltypeXexp = {};
    FdbOn_ERP_blocksecXtrltypeXexp{1,1,1} =  FdbOn_ERP_blockbegin_punish_unexp;

    %Over blockSection x TrialType
    %try avg over elements
    ERP_blockbegin_punish_all(sbj,:,:) = (FdbOn_ERP_blockbegin_punish_unexp(sbj,:,:) + FdbOn_ERP_blockbegin_punish_exp(sbj,:,:))/2;
    ERP_blockbegin_reward_all(sbj,:,:) = (FdbOn_ERP_blockbegin_reward_unexp(sbj,:,:) + FdbOn_ERP_blockbegin_reward_exp(sbj,:,:))/2;
    ERP_blockmid_punish_all(sbj,:,:) = (FdbOn_ERP_blockmid_punish_unexp(sbj,:,:) + FdbOn_ERP_blockmid_punish_exp(sbj,:,:))/2;
    ERP_blockmid_reward_all(sbj,:,:) = (FdbOn_ERP_blockmid_reward_unexp(sbj,:,:) + FdbOn_ERP_blockmid_reward_exp(sbj,:,:))/2;
    ERP_blockfin_punish_all(sbj,:,:) = (FdbOn_ERP_blockfin_punish_unexp(sbj,:,:) + FdbOn_ERP_blockfin_punish_exp(sbj,:,:))/2;
    ERP_blockfin_reward_all(sbj,:,:) = (FdbOn_ERP_blockfin_reward_unexp(sbj,:,:) + FdbOn_ERP_blockfin_reward_exp(sbj,:,:))/2;

    %blockSection x Feedback (expectation)
    FdbOn_ERP_blocksecXfdb = {};
    FdbOn_ERP_blockbegin_unexp_all(sbj,:,:) = (FdbOn_ERP_blockbegin_punish_unexp(sbj,:,:) + FdbOn_ERP_blockbegin_reward_unexp(sbj,:,:))/2;
    FdbOn_ERP_blockbegin_exp_all(sbj,:,:) = (FdbOn_ERP_blockbegin_punish_exp(sbj,:,:) + FdbOn_ERP_blockbegin_reward_exp(sbj,:,:))/2;
    FdbOn_ERP_blockmid_unexp_all(sbj,:,:) = (FdbOn_ERP_blockmid_punish_unexp(sbj,:,:) + FdbOn_ERP_blockmid_reward_unexp(sbj,:,:))/2;
    FdbOn_ERP_blockmid_exp_all(sbj,:,:) = (FdbOn_ERP_blockmid_punish_exp(sbj,:,:) + FdbOn_ERP_blockmid_reward_exp(sbj,:,:))/2;
    FdbOn_ERP_blockfin_unexp_all(sbj,:,:) = (FdbOn_ERP_blockfin_punish_unexp(sbj,:,:) + FdbOn_ERP_blockfin_reward_unexp(sbj,:,:))/2;
    FdbOn_ERP_blockfin_exp_all(sbj,:,:) = (FdbOn_ERP_blockfin_punish_exp(sbj,:,:) + FdbOn_ERP_blockfin_reward_exp(sbj,:,:))/2;


    %Time across blocks
    block1 = block_idx(1,:);
    block2 = block_idx(2,:);
    block3 = block_idx(3,:);
    block1 = find(block1 ~= 0);
    block2 = find(block2 ~= 0);
    block3 = find(block3 ~= 0);
    
    FdbOn_ERP_allbegin_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(1,:),Unexpected_all_trial_idx)),3)); 
    FdbOn_ERP_allbegin_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(1,:),Expected_all_trial_idx)),3)); 
    FdbOn_ERP_allmid_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(2,:),Unexpected_all_trial_idx)),3)); 
    FdbOn_ERP_allmid_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(2,:),Expected_all_trial_idx)),3)); 
    FdbOn_ERP_allfin_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(3,:),Unexpected_all_trial_idx)),3)); 
    FdbOn_ERP_allfin_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(3,:),Expected_all_trial_idx)),3));


    %Within blocks: 16-sections
    num_sec = 16;
    block_sec_idx = cell(num_sec,1);            % section_number --> trial indices  
    sec_len = ceil(block_len / num_sec);
    trl = 1;
    %get trial indices
    for b = 1:num_blocks
        for s = 1:num_sec
            %sec_pos = 1;
            start = (b-1)*block_len + (s-1)*sec_len + 1;
            fin = (b-1)*block_len + s*sec_len;
                for t = 1:sec_len
                    if trl <= length(all_idx)                 %may not have all sections of a block (ex. if end was problematic)
                        if (start <= all_idx(trl) & all_idx(trl) <= fin)
                            %block_sec_idx{s}(sec_pos) = all_idx(trl);
                            %sec_pos = sec_pos + 1;
                            block_sec_idx{s} = cat(2, block_sec_idx{s}, all_idx(trl));
                            trl = trl+1;
                        end
                    end
                end
        end
    end
   
    fprintf("Defining trial indices took %d seconds", toc);
    

    
    %FRN calc
    toi = [.228 .344];
    toi_idx = dsearchn(time',toi');                                     
    ERP_blocksecs = cell(num_sec,4);   %Section   {1}ERP_unexp   {2}ERP_exp  {3}Diff_Wave  {4}FRN
    for s = 1:num_sec
        ERP_blocksecs{s,1} = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_sec_idx{s},Unexpected_all_trial_idx)),3)); 
        ERP_blocksecs{s,2} = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_sec_idx{s},Expected_all_trial_idx)),3)); 
        %DiffWave: unexpected (neg) - expected(pos)
        ERP_blocksecs{s,3} = ERP_blocksecs{s,1} - ERP_blocksecs{s,2};
        %FRN: toi - baseline
        ERP_blocksecs{s,4} = mean(ERP_blocksecs{s,3}(choi_idx,toi_idx(1):toi_idx(2)) - mean(ERP_blocksecs{s,3}(choi_idx,baseline_idx(1):baseline_idx(2))));
    end
    %add to structures for all sbj!
    DiffWave_blocksecs{sbj} = ERP_blocksecs {:,3};   %all sections, DiffWave
    FRN_blocksec(sbj,:) = ERP_blocksecs(:,4)';         %all sections, FRN


    
    
    %%
    FdbOn_ERP_allbegin_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(1,:),Expected_all_trial_idx)),3)); 
    FdbOn_ERP_allmid_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(2,:),Unexpected_all_trial_idx)),3)); 
    FdbOn_ERP_allmid_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(2,:),Expected_all_trial_idx)),3)); 
    FdbOn_ERP_allfin_unexp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(3,:),Unexpected_all_trial_idx)),3)); 
    FdbOn_ERP_allfin_exp(sbj,:,:) = squeeze(mean(eeg_data(EEG_chan_indices,:,intersect(block_idx(3,:),Expected_all_trial_idx)),3));

    %FRN_pipeline
    %for i = 1:all_idx
    %FdbOn_ERP_alltrials_data(sbj,:,:i) = 

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
    FRN_Reward(sbj) = mean(DiffWave_Reward(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_Reward(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));

    DiffWave_Punish(sbj,:,:) = FdbOn_ERP_Punish_neg(sbj,:,:) - FdbOn_ERP_Punish_pos(sbj,:,:);
    FRN_Punish(sbj) = mean(DiffWave_Punish(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_Punish(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));


    %DiffWave : change over time! ----

    %Within Blocks: block section X fdb expectation
    DiffWave_Blockbegin(sbj,:,:) = FdbOn_ERP_blockbegin_unexp_all(sbj,:,:) - FdbOn_ERP_blockbegin_exp_all(sbj,:,:);
    FRN_Blockbegin(sbj) = mean(DiffWave_Blockbegin(sbj, choi_idx, toi_idx(1):toi_idx(2))) - mean(DiffWave_Blockbegin(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));

    DiffWave_Blockmid(sbj,:,:) = FdbOn_ERP_blockmid_unexp_all(sbj,:,:) - FdbOn_ERP_blockmid_exp_all(sbj,:,:);
    FRN_Blockmid(sbj) = mean(DiffWave_Blockmid(sbj, choi_idx, toi_idx(1):toi_idx(2))) - mean(DiffWave_Blockmid(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));

    DiffWave_Blockfin(sbj,:,:) = FdbOn_ERP_blockfin_unexp_all(sbj,:,:) - FdbOn_ERP_blockfin_exp_all(sbj,:,:);
    FRN_Blockfin(sbj) = mean(DiffWave_Blockfin(sbj, choi_idx, toi_idx(1):toi_idx(2)))- mean(DiffWave_Blockfin(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));


    %Across Blocks: block num X fdb expectation
    DiffWave_allbegin(sbj,:,:) = FdbOn_ERP_allbegin_unexp(sbj,:,:) - FdbOn_ERP_allbegin_exp(sbj,:,:);
    FRN_allbegin(sbj) = mean(DiffWave_allbegin(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_allbegin(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));

    DiffWave_allmid(sbj,:,:) = FdbOn_ERP_allmid_unexp(sbj,:,:) - FdbOn_ERP_allmid_exp(sbj,:,:);
    FRN_allmid(sbj) = mean(DiffWave_allmid(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_allmid(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));


    DiffWave_allfin(sbj,:,:) = FdbOn_ERP_allfin_unexp(sbj,:,:) - FdbOn_ERP_allfin_exp(sbj,:,:);
    FRN_allfin(sbj) = mean(DiffWave_allfin(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_allfin(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));


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
    
        %get eeg_data for trial_type x expectation
        %(looking at fdblocked only)
        eeg_data_fdbLocked_punish_unexp = eeg_data(choi_idx,:,Punish_unexpected_trial_idx); 
        eeg_data_fdbLocked_punish_exp = eeg_data(choi_idx,:,Punish_expected_trial_idx);
        eeg_data_fdbLocked_reward_unexp = eeg_data(choi_idx,:,Reward_unexpected_trial_idx);
        eeg_data_fdbLocked_reward_exp = eeg_data(choi_idx,:, Reward_expected_trial_idx);

        %get eeg_data for time x trial_type x expectation
        eeg_data_block_begin_punish_unexp = eeg_data(choi_idx,:,block_begin_punish_unexp_idx);
        eeg_data_block_begin_punish_exp = eeg_data(choi_idx,:,block_begin_punish_exp_idx);
        eeg_data_block_begin_reward_unexp = eeg_data(choi_idx,:,block_begin_reward_unexp_idx);
        eeg_data_block_begin_reward_exp = eeg_data(choi_idx,:,block_begin_reward_exp_idx);

        eeg_data_block_mid_punish_unexp = eeg_data(choi_idx,:,block_mid_punish_unexp_idx);
        eeg_data_block_mid_punish_exp = eeg_data(choi_idx,:,block_mid_punish_exp_idx);
        eeg_data_block_mid_reward_unexp = eeg_data(choi_idx,:,block_mid_reward_unexp_idx);
        eeg_data_block_mid_reward_exp = eeg_data(choi_idx,:,block_mid_reward_exp_idx);

        eeg_data_block_fin_punish_unexp = eeg_data(choi_idx,:,block_fin_punish_unexp_idx);
        eeg_data_block_fin_punish_exp = eeg_data(choi_idx,:,block_fin_punish_exp_idx);
        eeg_data_block_fin_reward_unexp = eeg_data(choi_idx,:,block_fin_reward_unexp_idx);
        eeg_data_block_fin_reward_exp = eeg_data(choi_idx,:,block_fin_reward_exp_idx);

        %Time: within blocks
        eeg_data_block_begin_all = eeg_data(choi_idx, :, block_begin);
        eeg_data_block_mid_all = eeg_data(choi_idx, :, block_mid);
        eeg_data_block_fin_all = eeg_data(choi_idx,:, block_fin);

        %Time: between blocks
        eeg_data_block1 = eeg_data(choi_idx,:,block1);
        eeg_data_block2 = eeg_data(choi_idx,:,block2);
        eeg_data_block3 = eeg_data(choi_idx,:,block3);


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

    %add for TrialType and Expectation
    fft_data_punish_unexp = fft(reshape(eeg_data_fdbLocked_punish_unexp,1,t*size(eeg_data_fdbLocked_punish_unexp,3)), nConv,2);   %add for trial type
    fft_data_punish_exp = fft(reshape(eeg_data_fdbLocked_punish_exp,1,t*size(eeg_data_fdbLocked_punish_exp,3)), nConv,2);
    fft_data_reward_unexp = fft(reshape(eeg_data_fdbLocked_reward_unexp,1,t*size(eeg_data_fdbLocked_reward_unexp,3)), nConv,2);
    fft_data_reward_exp = fft(reshape(eeg_data_fdbLocked_reward_exp,1,t*size(eeg_data_fdbLocked_reward_exp,3)), nConv,2);

    %fft for Time x TrialType x Expectation
    fft_data_begin_punish_unexp = fft(reshape(eeg_data_block_begin_punish_unexp,1,t*size(eeg_data_block_begin_punish_unexp,3)), nConv,2); 
    fft_data_begin_punish_exp = fft(reshape(eeg_data_block_begin_punish_exp,1,t*size(eeg_data_block_begin_punish_exp,3)), nConv,2); 
    fft_data_begin_reward_unexp = fft(reshape(eeg_data_block_begin_reward_unexp,1,t*size(eeg_data_block_begin_reward_unexp,3)), nConv,2); 
    fft_data_begin_reward_exp = fft(reshape(eeg_data_block_begin_reward_exp,1,t*size(eeg_data_block_begin_reward_exp,3)), nConv,2); 

    fft_data_mid_punish_unexp = fft(reshape(eeg_data_block_mid_punish_unexp,1,t*size(eeg_data_block_mid_punish_unexp,3)), nConv,2); 
    fft_data_mid_punish_exp = fft(reshape(eeg_data_block_mid_punish_exp,1,t*size(eeg_data_block_mid_punish_exp,3)), nConv,2); 
    fft_data_mid_reward_unexp = fft(reshape(eeg_data_block_mid_reward_unexp,1,t*size(eeg_data_block_mid_reward_unexp,3)), nConv,2); 
    fft_data_mid_reward_exp = fft(reshape(eeg_data_block_mid_reward_exp,1,t*size(eeg_data_block_mid_reward_exp,3)), nConv,2); 

    fft_data_fin_punish_unexp = fft(reshape(eeg_data_block_fin_punish_unexp,1,t*size(eeg_data_block_fin_punish_unexp,3)), nConv,2); 
    fft_data_fin_punish_exp = fft(reshape(eeg_data_block_fin_punish_exp,1,t*size(eeg_data_block_fin_punish_exp,3)), nConv,2); 
    fft_data_fin_reward_unexp = fft(reshape(eeg_data_block_fin_reward_unexp,1,t*size(eeg_data_block_fin_reward_unexp,3)), nConv,2); 
    fft_data_fin_reward_exp = fft(reshape(eeg_data_block_fin_reward_exp,1,t*size(eeg_data_block_fin_reward_exp,3)), nConv,2); 

    %fft for Time: block section 
    fft_data_fin_all = fft(reshape(eeg_data_block_fin_all,1,t*size(eeg_data_block_fin_all,3)), nConv,2); 
    fft_data_mid_all = fft(reshape(eeg_data_block_mid_all,1,t*size(eeg_data_block_mid_all,3)), nConv,2); 
    fft_data_begin_all = fft(reshape(eeg_data_block_begin_all,1,t*size(eeg_data_block_begin_all,3)), nConv,2); 

    %Time: between blocks
    fft_data_block1 = fft(reshape(eeg_data_block1,1,t*size(eeg_data_block1,3)), nConv,2); 
    fft_data_block2 = fft(reshape(eeg_data_block2,1,t*size(eeg_data_block2,3)), nConv,2); 
    fft_data_block3 = fft(reshape(eeg_data_block3,1,t*size(eeg_data_block3,3)), nConv,2); 





    bsl_win_TF = [-.4, -.2]; %baseline time window for TF analysis

    bsl_st_idx = dsearchn(time',bsl_win_TF(1));%idx of baseline window start
    bsl_end_idx = dsearchn(time',bsl_win_TF(1));%idx of baseline window end

    %Initiate matrix for TF power (for indiv sbj) 
    TF_Power_StimLocked = zeros(nFrex,npnts);
    TF_Power_FdbLocked = zeros(nFrex,npnts);

    %add matrices for Trial Type x Expectation
    TF_Power_Punish_unexp = zeros(nFrex,npnts);
    TF_Power_Punish_exp = zeros(nFrex,npnts);
    TF_Power_Reward_unexp = zeros(nFrex,npnts);
    TF_Power_Reward_exp = zeros(nFrex,npnts);

    %add matrices for Time x TrialType x Expectation ---
%     TF_Power_begin_punish_unexp;

    %add matrices for Time: within blocks 
    TF_Power_Begin = zeros(nFrex,npnts);
    TF_Power_Mid = zeros(nFrex,npnts);
    TF_Power_Fin = zeros(nFrex,npnts);

    %add matrices for Time: between blocks
    TF_Power_Block1 = zeros(nFrex,npnts);
    TF_Power_Block2 = zeros(nFrex,npnts);
    TF_Power_Block3 = zeros(nFrex,npnts);
        


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

        %TrialType x Expectation 
        as = ifft(fft_data_punish_unexp .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_punish_unexp = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_punish_unexp(bsl_st_idx:bsl_end_idx));
        TF_Power_Punish_unexp(fi,:) = 10*log10(tmp_TF_Power_punish_unexp./bsl);% dB conversion
        
        as = ifft(fft_data_punish_exp .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_punish_exp = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_punish_exp(bsl_st_idx:bsl_end_idx));
        TF_Power_Punish_exp(fi,:) = 10*log10(tmp_TF_Power_punish_exp./bsl);% dB conversion

        as = ifft(fft_data_reward_unexp .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_reward_unexp = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_reward_unexp(bsl_st_idx:bsl_end_idx));
        TF_Power_Reward_unexp(fi,:) = 10*log10(tmp_TF_Power_reward_unexp./bsl);% dB conversion

        as = ifft(fft_data_reward_exp .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_reward_exp = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_reward_exp(bsl_st_idx:bsl_end_idx));
        TF_Power_Reward_exp(fi,:) = 10*log10(tmp_TF_Power_reward_exp./bsl);% dB conversion

        %Time(within block)
        as = ifft(fft_data_begin_all .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_begin = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_begin(bsl_st_idx:bsl_end_idx));
        TF_Power_Begin(fi,:) = 10*log10(tmp_TF_Power_begin./bsl);% dB conversion
        
        as = ifft(fft_data_mid_all .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_mid = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_mid(bsl_st_idx:bsl_end_idx));
        TF_Power_Mid(fi,:) = 10*log10(tmp_TF_Power_mid./bsl);% dB conversion

        as = ifft(fft_data_fin_all .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_fin = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_fin(bsl_st_idx:bsl_end_idx));
        TF_Power_Fin(fi,:) = 10*log10(tmp_TF_Power_fin./bsl);% dB conversion


        %Time Across Blocks
        as = ifft(fft_data_block1 .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_block1 = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_block1(bsl_st_idx:bsl_end_idx));
        TF_Power_Block1(fi,:) = 10*log10(tmp_TF_Power_block1./bsl);% dB conversion
        
        as = ifft(fft_data_block2 .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_block2 = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_block2(bsl_st_idx:bsl_end_idx));
        TF_Power_Block2(fi,:) = 10*log10(tmp_TF_Power_block2./bsl);% dB conversion

        as = ifft(fft_data_block3 .* cmwX);
        as = as(:,halfwav:end-halfwav+1);
        sig = reshape(as,t,ntrials)';
        tmp_TF_Power_block3 = mean(abs(sig).^2,1);
        bsl = mean(tmp_TF_Power_block3(bsl_st_idx:bsl_end_idx));
        TF_Power_Block3(fi,:) = 10*log10(tmp_TF_Power_block3./bsl);% dB conversion



    end

    TF_Pow_StimLocked{sbj} = TF_Power_StimLocked;
    TF_Pow_FdbLocked{sbj} = TF_Power_FdbLocked;

    TF_Pow_Punish_unexp{sbj} = TF_Power_Punish_unexp;
    TF_Pow_Punish_exp{sbj} = TF_Power_Punish_exp;
    TF_Pow_Reward_unexp{sbj} = TF_Power_Reward_unexp;
    TF_Pow_Reward_exp{sbj} = TF_Power_Reward_exp;

    %Time: within blocks
    TF_Pow_Begin{sbj} = TF_Power_Begin;
    TF_Pow_Mid{sbj} = TF_Power_Mid;
    TF_Pow_Fin{sbj} = TF_Power_Fin;

    %Time: between blocks
    TF_Pow_Blk1{sbj} = TF_Power_Block1;
    TF_Pow_Blk2{sbj} = TF_Power_Block2;
    TF_Pow_Blk3{sbj} = TF_Power_Block3;

    %% FRN & Time

%     FdbOn_ERP_block_begin_all(sbj,:,:) = squeeze(mean(eeg_data_block_begin_all),3)); 
%     FdbOn_ERP_block_mid_all(sbj,:,:) = squeeze(mean(eeg_data_block_mid_all),3);
%     FdbOn_ERP_block_fin_all(sbj,:,:) = squeeze(mean(eeg_data_block_fin_all),3);

    %DiffWave_Punish(sbj,:,:) = FdbOn_ERP_Punish_neg(sbj,:,:) - FdbOn_ERP_Punish_pos(sbj,:,:);
    %FRN_Punish(sbj) = mean(DiffWave_Punish(sbj,choi_idx,toi_idx(1):toi_idx(2))) - mean(DiffWave_Punish(sbj,choi_idx,baseline_idx(1):baseline_idx(2)));

    

    

%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% plotting individual subjects
% 
%     if plot_indv
%         %plotting difference wave
%         fig=figure('Position', [80 80 800 400]);
%         tl = tiledlayout(1,2);
%         tl.TileSpacing = 'compact';
%         title(tl,['Sbj # ' num2str(sbj)])
% 
%         nexttile
%         plot(time,squeeze(FdbOn_ERP_pos(sbj,choi_idx,:)),'b');hold on;
%         plot(time,squeeze(FdbOn_ERP_neg(sbj,choi_idx,:)),'r');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title(' Pos Fdb ERP vs Neg Fdb ERP')
% 
%         nexttile
%         plot(time,squeeze(DiffWave(sbj,choi_idx,:)),'b');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title(' Difference Wave')
% 
%         mkdir(fullfile(figFolder,'Individual Sbj Plots'))
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave.jpg']));
%         close(fig)
% 
%         %plotting difference wave For Reward Trials 
%         fig=figure('Position', [80 80 800 400]);
%         tl = tiledlayout(1,2);
%         tl.TileSpacing = 'compact';
%         title(tl,['Sbj # ' num2str(sbj)])
% 
%         nexttile
%         plot(time,squeeze(FdbOn_ERP_Reward_pos(sbj,choi_idx,:)),'b');hold on;
%         plot(time,squeeze(FdbOn_ERP_Reward_neg(sbj,choi_idx,:)),'r');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title({' Pos Fdb ERP vs Neg Fdb ERP', 'Reward Trials'})
% 
%         nexttile
%         plot(time,squeeze(DiffWave_Reward(sbj,choi_idx,:)),'b');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title({'Difference Wave', 'Reward Trials'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave_Reward.jpg']));
%         close(fig)
% 
%         %plotting difference wave For Punish Trials 
%         fig=figure('Position', [80 80 800 400]);
%         tl = tiledlayout(1,2);
%         tl.TileSpacing = 'compact';
%         title(tl,['Sbj # ' num2str(sbj)])
% 
%         nexttile
%         plot(time,squeeze(FdbOn_ERP_Punish_pos(sbj,choi_idx,:)),'b');hold on;
%         plot(time,squeeze(FdbOn_ERP_Punish_neg(sbj,choi_idx,:)),'r');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         legend({'Pos Fdb','Neg Fdb'},'Location','northwest');
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title({' Pos Fdb ERP vs Neg Fdb ERP', 'Punish Trials'})
% 
%         nexttile
%         plot(time,squeeze(DiffWave_Punish(sbj,choi_idx,:)),'b');
%         xlim([-.2 .5])
%         ylabel('Amplitude [\mu V]')
%         xlabel('Time post Fdb [ms]')
%         x = [toi(1) toi(1) toi(2) toi(2)];
%         y = [min(ylim(gca)) min(ylim(gca))+.2 min(ylim(gca))+.2 min(ylim(gca))];
%         b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
%         b.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         X_zeroLine = xline(0,'k');
%         X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         Y_zeroLine = yline(0,'k');
%         Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         title({'Difference Wave', 'Punish Trials'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_DifferenceWave_Punish.jpg']));
%         close(fig)
% 
% 
%         %plot topoplots for Pos fdb Neg fdb and DiffWave
%         fig=figure('Position', [80 80 800 400],'Name',['Sbj' num2str(sbj)]);
%         tl = tiledlayout(1,3);
%         tl.TileSpacing = 'compact';
% 
%         nexttile
%         temp_data = squeeze(mean(FdbOn_ERP_neg(sbj,:,toi_idx(1):toi_idx(2)),3));
%         topoplot(temp_data',Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%         colormap('turbo')
%         c = colorbar;
%         c.Label.String = '\mu V';
%         title(['Neg Fdb'])
%         caxis([-5 5]);
% 
%         nexttile
%         temp_data = squeeze(mean(FdbOn_ERP_pos(sbj,:,toi_idx(1):toi_idx(2)),3));
%         topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%         colormap('turbo')
%         c = colorbar;
%         c.Label.String = '\mu V';
%         title(['Pos Fdb'])
%         caxis([-5 5]);
% 
%         nexttile
%         temp_data = squeeze(mean(DiffWave(sbj,:,toi_idx(1):toi_idx(2)),3));
%         topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%         colormap('turbo')
%         c = colorbar;
%         c.Label.String = '\mu V';
%         title(['Neg-Pos'])
%         caxis([-5 5]);
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_FRN_Topoplots.jpg']));
%         close(fig)
% 
%         %TF-plot for visualizeing FMT
%         fig=figure;
%         contourf(time,frex,TF_Power_FdbLocked,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Fdb Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_FdbLocked.jpg']));
%         close(fig)
% 
% 
%         fig=figure;
%         contourf(time,frex,TF_Power_StimLocked,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Stim Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Stimulus Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_StimLocked.jpg']));
%         close(fig)

%         %TrialType X Expectation 
%         fig=figure;
%         contourf(time,frex,TF_Power_Punish_exp,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Fdb Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_PunishTrls_ExpectedFdb.jpg']));
%         close(fig)

%         fig=figure;
%         contourf(time,frex,TF_Power_Punish_unexp,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Fdb Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_PunishTrls_UnexpectedFdb.jpg']));
%         close(fig)

%         fig=figure;
%         contourf(time,frex,TF_Power_Reward_exp,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Fdb Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_RewardTrls_ExpectedFdb.jpg']));
%         close(fig)
%        
%         fig=figure;
%         contourf(time,frex,TF_Power_Reward_unexp,40,'linecolor','none');
%         set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%         colormap(turbo)
%         ylabel('Freq [Hz]')
%         xlabel('Time from Fdb Onset')
%         xlim([-.5 2])
%         caxis([-5 5])
%         a=colorbar;
%         ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%         a.Label.Position(1) = 3;
%         a.Ticks = [-5,0,5];
%         title({['Sbj #' num2str(sbj) ' TF-Power'], 'Feedback Locked - FCz'})
% 
%         saveas(fig,fullfile(figFolder,'Individual Sbj Plots',['Sbj_' filenames(sbj).name(1:end-4) '_TF-Power_FCz_RewardTrls_UnexpectedFdb.jpg']));
%         close(fig)
%     end
end

%% Theta - Time within block
%gather FT data from structure
    theta_frq = 4:8;
    fmt_start_idx = dsearchn(time',0.200);    %(Eisma et al, 2021): uses time window of 200ms to 450ms after stimulus/reward
    fmt_end_idx = dsearchn(time',0.400);      %(Verbeke and Verguts, 2019):time window 250-500ms
                                              % (Cavanagh et al, 2010): 200-400ms


    %sbj X Hz X time point in trial
    FMT_data_punish_unexp = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Punish_unexp{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_punish_unexp(i,:,:) = FMT_data;
    end

    FMT_data_punish_exp = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Punish_exp{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_punish_exp(i,:,:) = FMT_data;
    end

    FMT_data_reward_unexp = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Reward_unexp{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_reward_unexp(i,:,:) = FMT_data;
    end

    FMT_data_reward_exp = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Reward_exp{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_reward_exp(i,:,:) = FMT_data;
    end



    %FMT change over TIME: within block
    %get data from indices
    FMT_data_begin = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Begin{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_begin(i,:,:) = FMT_data;
    end
    FMT_data_mid = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Mid{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_mid(i,:,:) = FMT_data;
    end
    FMT_data_fin = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Fin{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_fin(i,:,:) = FMT_data;
    end

    %% Time Across Blocks
    FMT_data_block1 = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Blk1{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_block1(i,:,:) = FMT_data;
    end
    FMT_data_block2 = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Blk2{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_block2(i,:,:) = FMT_data;
    end
    FMT_data_block3 = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
    for i = 1:num_sbj
        sbj_cellarr = TF_Pow_Blk3{i};
        FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
        FMT_data_block3(i,:,:) = FMT_data;
    end


%% Plotting Group Average
if plot_grandAvg
%     %plotting difference wave
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
        legend({'Exp Fdb','Unexp Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Expected Fdb ERP vs Unexpected Fdb ERP')
    
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
        title(tl,'Reward Trials')
    
        nexttile
        plot(time,squeeze(mean(FdbOn_ERP_Reward_pos(:,choi_idx,:),1)),'b');hold on;
        plot(time,squeeze(mean(FdbOn_ERP_Reward_neg(:,choi_idx,:),1)),'r');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        legend({'Expected Fdb','Unexpected Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({'Expected Fdb ERP vs Unexpected Fdb ERP'})
    
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
        title({'Difference Wave'})
    
        saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave_Reward.jpg'));
        close(fig)
    
        %plotting difference Punish Trials
        fig=figure('Position', [80 80 800 400]);
        tl = tiledlayout(1,2);
        tl.TileSpacing = 'compact';
        title(tl,'Punish Trials')
    
        nexttile
        plot(time,squeeze(mean(FdbOn_ERP_Punish_pos(:,choi_idx,:),1)),'b');hold on;
        plot(time,squeeze(mean(FdbOn_ERP_Punish_neg(:,choi_idx,:),1)),'r');
        xlim([-.2 .5])
        ylabel('Amplitude [\mu V]')
        xlabel('Time post Fdb [ms]')
        legend({'Expected Fdb','Unexpected Fdb'},'Location','northwest');
        x = [toi(1) toi(1) toi(2) toi(2)];
        y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
        b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
        b.Annotation.LegendInformation.IconDisplayStyle = 'off';
        X_zeroLine = xline(0,'k');
        X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Y_zeroLine = yline(0,'k');
        Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title({'Expected Fdb ERP vs Unexpected Fdb ERP'})
    
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
        title({'Difference Wave'})
    
        saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave_Punish.jpg'));
        close(fig)

      %plotting Diff Wave of ERPs
      DiffWave_TrlTypeXExp = DiffWave_Punish - DiffWave_Reward; 
      fig=figure('Position', [80 80 800 400]);
      tl = tiledlayout(1,2);
      tl.TileSpacing = 'compact';
      title(tl,'Grand Average')

      nexttile
      plot(time,squeeze(mean(DiffWave_Punish(:,choi_idx,:),1)),'b');hold on;
      plot(time,squeeze(mean(DiffWave_Reward(:,choi_idx,:),1)),'r');
      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Punish','Reward'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title({'Expectation Difference Wave (U-E)', 'Trial Type'})

      nexttile
      plot(time,squeeze(mean(DiffWave_TrlTypeXExp(:,choi_idx,:),1)),'b');
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
      title({'Difference Wave (P-R)', 'Trial Type'})

      saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_DifferenceWave_TrlType.jpg'));
      close(fig)

      %Block Section (time) Difference Wave (Unexp - Expected)
      fig=figure('Position', [80 80 800 400]);
      tl = tiledlayout(2,2);
      tl.TileSpacing = 'compact';
      title(tl,'Block Section Grand Average')

      nexttile
      plot(time,squeeze(mean(DiffWave_Blockbegin(:,choi_idx,:),1)),'r');hold on;
      plot(time,squeeze(mean(DiffWave_Blockmid(:,choi_idx,:),1)),'g');hold on;
      plot(time,squeeze(mean(DiffWave_Blockfin(:,choi_idx,:),1)),'b');

      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Begin','Mid', 'Fin'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title('RewP Difference Waves Across Block Sections')

      nexttile
      plot(time,squeeze(mean(DiffWave_Blockbegin(:,choi_idx,:),1)),'r');
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
      title('Section 1')

      nexttile
      plot(time,squeeze(mean(DiffWave_Blockmid(:,choi_idx,:),1)),'g');
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
      title('Section 2')

      nexttile
      plot(time,squeeze(mean(DiffWave_Blockfin(:,choi_idx,:),1)),'b');
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
      title('Section 3')

      mkdir(fullfile(figFolder,'Grand Average'))
      saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_BlockSect.jpg'));
      close(fig)


      %Time: Across Blocks RewP (Unexp - Expected)
      fig=figure('Position', [80 80 800 400]);
      tl = tiledlayout(2, 2);
      tl.TileSpacing = 'compact';
      title(tl,'RewP Across Blocks Grand Average')

      nexttile
      plot(time,squeeze(mean(DiffWave_allbegin(:,choi_idx,:),1)),'r');hold on;
      plot(time,squeeze(mean(DiffWave_allmid(:,choi_idx,:),1)),'g');hold on;
      plot(time,squeeze(mean(DiffWave_allfin(:,choi_idx,:),1)),'b');
      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Begin','Mid', 'Fin'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title('RewP Difference Waves Across Blocks')

      nexttile
      plot(time,squeeze(mean(DiffWave_allbegin(:,choi_idx,:),1)),'r');
      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Begin'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title('RewP During Block 1')

      nexttile
      plot(time,squeeze(mean(DiffWave_allmid(:,choi_idx,:),1)),'g');
      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Begin','Mid', 'Fin'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title('RewP During Block 2')

      nexttile
      plot(time,squeeze(mean(DiffWave_allfin(:,choi_idx,:),1)),'b');
      xlim([-.2 .5])
      ylabel('Amplitude [\mu V]')
      xlabel('Time post Fdb [ms]')
      legend({'Begin','Mid', 'Fin'},'Location','northwest');
      x = [toi(1) toi(1) toi(2) toi(2)];
      y = [min(ylim(gca)) min(ylim(gca))+.25 min(ylim(gca))+.25 min(ylim(gca))];
      b = patch(x,y,'r','FaceAlpha',.25,'LineStyle','none');
      b.Annotation.LegendInformation.IconDisplayStyle = 'off';
      X_zeroLine = xline(0,'k');
      X_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      Y_zeroLine = yline(0,'k');
      Y_zeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
      title('RewP During Block 3')

      mkdir(fullfile(figFolder,'Grand Average'))
      saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_AcrossBlocks.jpg'));
      close(fig)



      %%plot TOPOPLOTS for Pos fdb Neg fdb and DiffWave
%     fig=figure('Position', [80 80 800 400],'Name','Grand Average');
%     tl = tiledlayout(1,3);
%     tl.TileSpacing = 'compact';
% 
%     nexttile
%     temp_data = squeeze(mean(mean(FdbOn_ERP_neg(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data',Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Neg Fdb'])
%     caxis([-5 5]);
% 
%     nexttile
%     temp_data = squeeze(mean(mean(FdbOn_ERP_pos(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Pos Fdb'])
%     caxis([-5 5]);
% 
%     nexttile
%     temp_data = squeeze(mean(mean(DiffWave(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Neg-Pos'])
%     caxis([-5 5]);
% 
%     saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_Topoplots.jpg'));
%     close(fig)

%     %TrialType X Expectation
%     fig=figure('Position', [80 80 800 400],'Name','Grand Average');
%     tl = tiledlayout(1,3);
%     tl.TileSpacing = 'compact';
% 
%     nexttile
%     temp_data = squeeze(mean(mean(FdbOn_ERP_neg(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data',Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Neg Fdb'])
%     caxis([-5 5]);
% 
%     nexttile
%     temp_data = squeeze(mean(mean(FdbOn_ERP_pos(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Pos Fdb'])
%     caxis([-5 5]);
% 
%     nexttile
%     temp_data = squeeze(mean(mean(DiffWave(:,:,toi_idx(1):toi_idx(2)),3),1));
%     topoplot(temp_data,Chan_locs_acticap_vWMExp2a,'plotrad',.53,'electrodes','on','style','map');
%     colormap('turbo')
%     c = colorbar;
%     c.Label.String = '\mu V';
%     title(['Neg-Pos'])
%     caxis([-5 5]);
% 
%     saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_FRN_Topoplots.jpg'));
%     close(fig)

   
     
      %% TF-plot for visualising FMT POWER
%     fig=figure;
%     contourf(time,frex,squeeze(mean(cat(3,TF_Pow_FdbLocked{:}),3)),40,'linecolor','none');
%     set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%     colormap(turbo)
%     ylabel('Freq [Hz]')
%     xlabel('Time from Fdb Onset')
%     xlim([-.5 2])
%     caxis([-5 5])
%     a=colorbar;
%     ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%     a.Label.Position(1) = 3;
%     a.Ticks = [-5,0,5];
%     title({'Grand Average TF-Power', 'Feedback Locked - FCz'})
% 
%     saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_FdbLocked.jpg'));
%     close(fig)
% 
% 
%     fig=figure;
%     contourf(time,frex,squeeze(mean(cat(3,TF_Pow_FdbLocked{:}),3)),40,'linecolor','none');
%     set(gca,'ytick',[4, 8, 13, 22, 30, 55],'yscale','log','YMinorTick','off');
%     colormap(turbo)
%     ylabel('Freq [Hz]')
%     xlabel('Time from Stim Onset')
%     xlim([-.5 2])
%     caxis([-5 5])
%     a=colorbar;
%     ylabel(a,'dB change from bsl','FontSize',10,'Rotation',270);
%     a.Label.Position(1) = 3;
%     a.Ticks = [-5,0,5];
%     title({'Grand Average TF-Power', 'Stimulus Locked - FCz'})
% 
%     saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_StimLocked.jpg'));
%     close(fig)

      %punish unexpected
      fig=figure;
      tl = tiledlayout(1,3);
      tl.TileSpacing = 'compact';
      contourf(time,frex,squeeze(mean(cat(3,TF_Pow_Punish_unexp{:}),3)),40,'linecolor','none');
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
      title({'Punish Trials Unexpected Fdb'})

%       saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_PunishTrls_UnexpectedFdb.jpg'));
%       close(fig)

      %punish expected
      nexttile;
      contourf(time,frex,squeeze(mean(cat(3,TF_Pow_Punish_exp{:}),3)),40,'linecolor','none');
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
      title({'Punish Trials Expected Fdb'})

%       saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_PunishTrls_ExpectedFdb.jpg'));
%       close(fig)

      %reward unexpected
%       fig=figure;
      nexttile;
      contourf(time,frex,squeeze(mean(cat(3,TF_Pow_Reward_unexp{:}),3)),40,'linecolor','none');
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
      title({'Reward Trials Unexpected Fdb'})

%       saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_RewardTrls_UnexpectedFdb.jpg'));
%       close(fig)

      %reward expected
%       fig=figure;
      nexttile;
      contourf(time,frex,squeeze(mean(cat(3,TF_Pow_Reward_exp{:}),3)),40,'linecolor','none');
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
      title({'Grand Average TF-Power', 'Reward Trials Expected Fdb - FCz'})

      saveas(fig,fullfile(figFolder,'Grand Average','GrandAvg_TF-Power_FCz_TrialCond_FdbExp.jpg'));
      close(fig)

      %plotting changes across time
      %RewP amplitude

end

%% Stats?

%RewP: Unexpected vs Expected fdb
%note: not truly a RewP/FRN bc use total ERP, not difference waves
%using sbj means:
FRN_Unexp = mean(FdbOn_ERP_neg(:,choi_idx,toi_idx(1):toi_idx(2)), [2 3]); 
FRN_Exp = mean(FdbOn_ERP_pos(:,choi_idx,toi_idx(1):toi_idx(2)), [2 3]);

[h_fdbExp, p_fdbExp, ci_fdbExp, stats_fdbExp] = ttest(FRN_Unexp, FRN_Exp);

%using raw data
FRN_Unexp_r = reshape(FdbOn_ERP_neg(:,choi_idx,toi_idx(1):toi_idx(2)), [], 1);
FRN_Exp_r = reshape(FdbOn_ERP_pos(:,choi_idx,toi_idx(1):toi_idx(2)), [], 1);
[h_fdbExp_r, p_fdbExp_r, ci_fdbExp_r, stats_fdbExp_r] = ttest(FRN_Unexp, FRN_Exp);

%or: FRN --> diff from zero?
[h_fdbExp_0, p_fdbExp_0, ci_fdbExp_0, stats_fdbExp_0] = ttest(FRN);



%FRN Expectation difference: Punish vs Reward

    %descriptive
    FRN_Punish_mean = mean(FRN_Punish,"omitnan");
    FRN_Reward_mean = mean(FRN_Reward,"omitnan");
    FRN_Punish_sig = std(FRN_Punish, "omitnan");
    FRN_Reward_sig = std(FRN_Reward,"omitnan");
    
    %significance
    significance_results = struct();
    %h, p, ci, stats
    %[h, p, ci, stats] = ttest2(FRN_Punish, FRN_Reward);
    [h_TrlType,p_trltype,ci_trltype,stats_trltype] = ttest(FRN_Punish, FRN_Reward);   %paired-sample t-test. h=boolean for hypothesis, p=p-value, ci=confidence-interval, stats=pooled-variance...
    %Ttest_result_PvR_expectation = [h_TrlType,p_trltype,ci_trltype,stats_trltype];
    significance_results = ("");

%RewP TIME: across block sections
    %descriptive
    FRN_Blockbegin_mean = mean(FRN_Blockbegin,"omitnan");
    FRN_Blockmid_mean = mean(FRN_Blockmid,"omitnan");
    FRN_Blockfin_mean = mean(FRN_Blockfin,"omitnan");
    FRN_Blockbegin_sig = std(FRN_Blockbegin, "omitnan");
    FRN_Blockmid_sig = std(FRN_Blockmid, "omitnan");
    FRN_Blockfin_sig = std(FRN_Blockfin, "omitnan");

    %significance
    [p_blcksecXtrlt, tbl_blcksecXtrlt, stats_blcksecXtrlt] = anova1([FRN_Blockbegin; FRN_Blockmid; FRN_Blockfin]');
    [h, p, ci, stats] = ttest(FRN_Blockbegin, FRN_Blockmid);

%RewP TIME: across blocks: 3blocks
    %descriptive
    FRN_acrosblk_stats = {};                                         %mean      %sd         
    FRN_allbegin_mean = mean(FRN_allbegin,"omitnan");        %begin
    FRN_allmid_mean = mean(FRN_allmid,"omitnan");            %mid
    FRN_allfin_mean = mean(FRN_allfin,"omitnan");            %mid
    FRN_allbegin_sig = std(FRN_allbegin, "omitnan");
    FRN_allmid_sig = std(FRN_allmid, "omitnan");
    FRN_allfin_sig = std(FRN_allfin, "omitnan");

    %significnce
    [p_acrosblk, tbl_acrosblk, stats_acrosblk] = anova1([FRN_allbegin; FRN_allmid; FRN_allfin]');
    FRN_Blockbegin_Blockmid_stats = {};
    [h_ab, p_ab, ci_ab, stats_ab] = ttest(FRN_allmid, FRN_Blockmid);

%RewP TIME: across blocks: custom (1/5)
    FRN_blocksec = cell2mat(FRN_blocksec);
    %significance
    [p_blocksecs, tbl_blocksecs, stats_blocksecs] = anova1(FRN_blocksec);

    %plot
    figure
    boxchart(FRN_blocksec, "notch", "on");    %notches capture 95% CI for median
    title("FRN across block sections");
    xlabel("Block Sections");
    ylabel("Amplitude (µV)");

    %-- try fitting linear model

%     Tbl_FRN_blocksec = table("Size", [num_sbj*num_sec, 3],'VariableNames', {'sbj', 'block section', 'RewP amplitude'});
%     Tbl_FRN_blocksec(:,3) = reshape(FRN_blocksec', [], 1);
%     Tbl_FRN_blocksec(:,2) = repmat((1:num_sec)', num_sbj,1);
%     Tbl_FRN_blocksec(:,1) = repelem((1:num_sbj)', num_sec,1);

    sbj_id = reshape(FRN_blocksec', [], 1);
    block_section = repmat((1:num_sec)', num_sbj,1);
    RewP_amplitude = repelem((1:num_sbj)', num_sec,1);
    Tbl_FRN_blocksec = table(RewP_amplitude, block_section, sbj_id, 'VariableNames', {'sbj', 'block section', 'RewP amplitude'});
   
    model_FRN_blocksec = fitlm(Tbl_FRN_blocksec);

    %try repeated measures anova
    %rm_FRN_blocksec = fitrm(Tbl_FRN_blocksec, "RewP amplitude ~ block section");  %fit repeated measures model

    %change format: sbj ... blocks
    Tbl_FRN_blocksec_ranova = array2table(FRN_blocksec);
    Tbl_FRN_blocksec_ranova.("sbj") = (1:num_sbj)';
    
    rm_FRN_blocksec = fitrm(Tbl_FRN_blocksec_ranova, "FRN_blocksec1-FRN_blocksec16 ~ sbj");  %fit repeated measures model
    
    ranovatbl_FRN_blocksec = ranova(rm_FRN_blocksec);

%---FMT---

%(1): Feedback vs Stimulus-locked (conirm soundness of data)
    %sbj X Hz X time point in trial
    FMT_data_FdbLocked = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
        for i = 1:num_sbj
            sbj_cellarr = TF_Pow_FdbLocked{i};
            FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
            FMT_data_FdbLocked(i,:,:) = FMT_data;
        end
    FMT_data_StimLocked = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1);
        for i = 1:num_sbj
            sbj_cellarr = TF_Pow_StimLocked{i};
            FMT_data = sbj_cellarr(theta_frq, fmt_start_idx:fmt_end_idx);
            FMT_data_StimLocked(i,:,:) = FMT_data;
        end
    FMT_FdbLocked_raw = reshape(FMT_data_FdbLocked,[],1);
    FMT_StimLocked_raw = reshape(FMT_data_StimLocked,[],1);
    [h_FMT_fdbVStim, p_FMT_fdbVStim, ci_FMT_fdbVStim, stats_FMT_fdbVStim] = ttest(FMT_FdbLocked_raw, FMT_StimLocked_raw);
    
    CohenD_FMT_fdbVStim = meanEffectSize(reshape(FMT_data_FdbLocked,[],1), reshape(FMT_data_StimLocked,[],1));
    
    %by subj means
    [h_FMT_fdbVStim2, p_FMT_fdbVStim2, ci_FMT_fdbVStim2, stats_FMT_fdbVStim2] = ttest(mean(FMT_data_FdbLocked,[2 3]), mean(FMT_data_StimLocked,[2 3]));


   
%Trial Condition X Fdb Exp
%descriptive
FMT_power_punish_unexp = mean(FMT_data_punish_unexp,"all");
FMT_power_punish_exp = mean(FMT_data_punish_exp,"all");
FMT_power_reward_unexp = mean(FMT_data_reward_unexp,"all");
FMT_power_reward_exp = mean(FMT_data_reward_exp,"all");

%reshape       %unexp    %exp
%  punish
%  reward 

FMT_power_cond_exp = [];  
reps = num_sbj * 5 * 201;  %sbj x freq x time points
FMT_power_cond_exp(:,1) = cat(1, reshape(FMT_data_punish_unexp,[],1), reshape(FMT_data_reward_unexp,[],1));   
FMT_power_cond_exp(:,2) = cat(1,reshape(FMT_data_punish_exp,[],1), reshape(FMT_data_reward_exp,[],1));

[p_FMT_condFdb, tbl_FMT_condFdb, stats_FMT_condFdb] = anova2(FMT_power_cond_exp, reps);
[c_FMT_condFdb, m_FMT_condFdb] = multcompare(stats_FMT_condFdb);

%sbj means
FMT_power_cond_exp_sbjm = [];
FMT_power_cond_exp_sbjm(:,1) = cat(1, mean(FMT_data_punish_unexp, [2 3]), mean(FMT_data_reward_unexp,[2 3]));
FMT_power_cond_exp_sbjm(:,2) = cat(1, mean(FMT_data_punish_exp, [2 3]), mean(FMT_data_reward_exp, [2 3])); 

[p_FMT_condFdb_sbj, tbl_FMT_condFdb_sbj, stats_FMT_condFdb_sbj] = anova2(FMT_power_cond_exp_sbjm, num_sbj);



%repeated measures ANOVA (ranova)
FMT_punish_unexp_grped = reshape(FMT_data_punish_unexp,[],1);
FMT_punish_exp_grped = reshape(FMT_data_reward_unexp,[],1);
FMT_reward_unexp_grped = reshape(FMT_data_punish_exp,[],1);
FMT_reward_exp_grped = reshape(FMT_data_reward_exp,[],1);
FMT_table = table(FMT_punish_unexp_grped, FMT_punish_exp_grped, FMT_reward_unexp_grped, FMT_reward_exp_grped, 'VariableNames',{'Punish Unexp','Punish Exp', 'Reward Unexp', 'Reward Exp'});

rm = fitrm(FMT_table,'FMT_punish_unexp_grped-FMT_reward_exp_grped ~ 1');

    %difference
    FMT_punish_diff = FMT_data_punish_unexp - FMT_data_punish_exp;
    FMT_reward_diff = FMT_data_reward_unexp - FMT_data_reward_exp;

    FMT_punish_diff_mean = squeeze(mean(FMT_punish_diff, "all"));
    FMT_reward_diff_mean = squeeze(mean(FMT_reward_diff, "all"));

    %sbj means
    FMT_punish_diff_sbjmean = mean(FMT_punish_diff, [2 3]);
    FMT_reward_diff_sbjmean = mean(FMT_reward_diff, [2 3]);

    %descriptive stats: sd
    FMT_sd_punish_diff_sbjmean = std(FMT_punish_diff_sbjmean);
    FMT_sd_reward_diff_sbjmean = std(FMT_reward_diff_sbjmean);

    %Difference scores (visual check)
    FMT_diff = FMT_punish_diff_sbjmean - FMT_reward_diff_sbjmean;

%significance
[h_fmt, p_fmt, ci_fmt, stats_fmt] = ttest(FMT_punish_diff_sbjmean, FMT_reward_diff_sbjmean);




%FMT Power over Time: within block
%FMT_data_fin = zeros(num_sbj, 5, fmt_end_idx-fmt_start_idx + 1); %8 x 5 x _
    %group by subject
    FMT_fin_sbj = mean(FMT_data_fin, [2 3]);
    FMT_mid_sbj = mean(FMT_data_mid, [2 3]);
    FMT_begin_sbj = mean(FMT_data_mid, [2 3]);

    %descriptive
    FMT_fin_mean = mean(FMT_data_fin,"all");
    FMT_mid_mean = mean(FMT_data_mid,"all");
    FMT_begin_mean = mean(FMT_data_begin,"all");

    FMT_fin_sd = std(FMT_fin_sbj);
    FMT_mid_sd = std(FMT_mid_sbj);
    FMT_begin_sd = std(FMT_begin_sbj);

    %ANOVA
    FMT_time_sbj = cat(2, FMT_fin_sbj,FMT_mid_sbj,FMT_begin_sbj);  %col for diff time conditions
    [p_FMTblockw, table_FMTblockw, stats_FMT_blockw] = anova1(FMT_time_sbj);

    %compute with raw data 
    FMT_fin_raw = reshape(FMT_data_fin,[],1);
    reps2 = size(FMT_fin_raw, 1);S
    FMT_block3rds = cat(2, FMT_fin_raw, reshape(FMT_data_mid,[],1), reshape(FMT_data_fin, [],1));
    
    [p_FMTblockw_r, table_FMTblockw_r, stats_FMT_blockw_r] = anova1(FMT_block3rds);


   %FMT Power over Time: across blocks

   %ANOVA
   FMT_block1_sbj = mean(FMT_data_block1, [2 3]);
   FMT_block2_sbj = mean(FMT_data_block2, [2 3]);
   FMT_block3_sbj = mean(FMT_data_block3, [2 3]);
   
   FMT_time2_sbj = cat(2, FMT_block1_sbj,FMT_block2_sbj,FMT_block3_sbj);  %col for diff time conditions
   [p_FMTblockz, table_FMTblockz, stats_FMTblockz] = anova1(FMT_time2_sbj);

   %v2: consider all data (sbj x Hz x trials), without averaging across sbj
   FMT_block1_sbj_2 = reshape(FMT_data_block1,[],1);
   FMT_block2_sbj_2 = reshape(FMT_data_block2, [], 1);
   FMT_block3_sbj_2 = reshape(FMT_data_block3, [], 1);
   
   FMT_time2_sbj_2 = cat(2, FMT_block1_sbj_2,FMT_block2_sbj_2,FMT_block3_sbj_2);  %col for diff time conditions
   [p_FMTblockz2, table_FMTblockz2, stats_FMTblockz2] = anova1(FMT_time2_sbj_2);

   [h_FMT_timez3, p_FMTblockz3, ci_FMTBlockz3, statsFMTBlockz3] = ttest(FMT_time2_sbj_2(:,1), FMT_time2_sbj_2(:,3));

        %raw data --> large df --> more likely to find p < alpha.
        %corrections?