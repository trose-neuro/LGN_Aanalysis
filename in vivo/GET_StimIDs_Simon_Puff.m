function ids = GET_StimIDs_Simon_Puff( AuxData, aux_samplingrate, level )


%  channels
Frames = AuxData(3,:);
Stims = AuxData(12,:);

frame_times = get_frame_times(Frames);
frame_times_level = frame_times(1:level:end);


StimOn = Stims(frame_times_level)>0.8;
stim_onsets  = find(diff(StimOn==1)>0);
stim_offsets = find(diff(StimOn==1)<0);

% chop his up into stim epochs - he used bursts of stim. Let's only detect
% the first pulse of burst

burst_on  = find(diff(stim_onsets)>20) + 1;
burst_off = find(diff(stim_offsets)>20);

burst_on  = [1 burst_on];
burst_off = [burst_off size(stim_offsets,2)];

trial_onsets  = stim_onsets(burst_on);
trial_offsets = stim_onsets(burst_off);

% TBD: 'll implement within burst stims some time later.

%debug figure
% figure;
% plot(StimOn); ylim([-.5 1.5]);
% vline(stim_onsets(burst_on),'r')
% vline(stim_offsets(burst_off),'k')


ids.StimBounds          = [trial_onsets' trial_offsets'];
ids.StimBoundsComment   = '[Airpuff onset Airpuff offset]';

% extract a few further stimparameters from aux_data
ids.ITI_duration  = median(diff(frame_times_level(stim_onsets(burst_on)))) / aux_samplingrate;
ids.stim_duration = median(diff(frame_times_level([trial_onsets' trial_offsets']),[],2)) / aux_samplingrate;
