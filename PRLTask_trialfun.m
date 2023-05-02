function trl = PRLTask_trialfun(cfg)
%last edited 4/18/2023


hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
event = event(strcmp('Stimulus', {event.type}));
if isempty(event)
    event = ft_read_event(cfg.dataset);
    event = event(strcmp('Stim', {event.type}));
end
values      = {event.value};
samples     = {event.sample};


% Trial starts and ends
t_starts    = find(strcmp(values,'S101'));
t_ends      = find(strcmp(values,'S102'));

% Cut the data into trials
for i = 1:length(t_starts)
    trials{i}.values    = values(t_starts(i):t_ends(i));
    trials{i}.samples   = samples(t_starts(i):t_ends(i));
end

% Identify trial and sample information
for i = 1:length(trials)
    trlnum(i)        = i;

    % trial type
    trial_type_rew  = 21;%reward trial marker
    trial_type_pun  = 22;%punishment trial marker
    trlType(i)      = find(ismember({'S 21','S 22'},trials{i}.values));
    
    % response type
    resp_opt        = 61;%optimal response marker
    resp_subopt     = 62;%suboptimal response marker
    resp_none       = 63;%no response 
    respType(i)     = find(ismember({'S 61','S 62','S 63'},trials{i}.values));

    % stim on 
    stim_on         = 32;
    StimOnSample(i)  = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 32'})));


    % feedback on
    fdb_on          = 36;
    fdbOnSample(i)  = cell2mat(trials{i}.samples(strcmp(trials{i}.values,{'S 36'})));

    if respType(i) ~= 3
        % feedback type
        fdb_pos         = 71; % positive feedback
        fdb_neg         = 72; % negative feedback
        fdb_neu         = 73; % neutral feedback
        fdbType(i)      = find(ismember({'S 71','S 72','S 73'},trials{i}.values));
    else
        fdbType(i)      = 4; % if the subject did not respond in this trial, then assign feedback type as 4. 
    end
end

% lock the data with respect to stim onset
prestim    = 1;
poststim   = 2;

pretrig    = -round(prestim * hdr.Fs);
posttrig   = round(poststim * hdr.Fs);

for i = 1:length(trials)
    begsample(i) = StimOnSample(i) + pretrig;
    endsample(i) = fdbOnSample(i) + posttrig;
    offset(i)    = pretrig;
    offset2(i)   = fdbOnSample(i) - begsample(i); %offset 2 gives the amount of samples between epoch start and feedback onset. This can be used to redefine trials relative to feedback onset.
end

trl = [begsample' endsample' offset' trlnum' trlType' respType' fdbType' StimOnSample' fdbOnSample' offset2'];

end

