function ids = GET_StimIDs_Chirp( AuxData, aux_samplingrate, level )
% TR2020

%  channels bscope2)
Frames  = AuxData(4,:);
Stims   = AuxData(8,:);
eye1    = AuxData(17,:);
eye2    = AuxData(18,:);

frame_times         = get_frame_times(Frames);
frame_times_level   = frame_times(1:level:end);

minsample_delta = 100; %minimum stim duration in level frames

StimOn = Stims(frame_times_level)>0.8;
Eye1On = eye1(frame_times_level)*-1+max(eye1)>0.8;
Eye2On = eye2(frame_times_level)>0.8;
Eye2On(end) = 1;
bino = Eye1On+Eye2On;

Eye1On_only = Eye1On - bino + 1;
Eye2On_only = Eye2On - bino + 1;
Eye2On_only(end) = 0;
Eye1On_only(end) = 0;

chirp_onsets_temp  = find(diff(StimOn==1)>0);
chirp_offsets_temp = find(diff(StimOn==1)<0);

bino_onsets_temp  = [1 find(diff(bino==1)<0)];
bino_offsets_temp  = find(diff(bino==1)>0);
 
bino_onsets  = bino_onsets_temp(find(bino_offsets_temp - bino_onsets_temp > minsample_delta));
bino_offsets = bino_offsets_temp(find(bino_offsets_temp - bino_onsets_temp > minsample_delta));

% generate cleaned bino binary
bino_clean = logical(zeros(size(bino)));

for i = 1:size(bino_onsets,2)
    bino_clean(bino_onsets(i):bino_offsets(i)) = 1;
end


% extract chirp stim on and offsets

chirp_on  = find(diff(chirp_onsets_temp)>minsample_delta) + 1;
chirp_off = find(diff(chirp_offsets_temp)>minsample_delta);
chirp_on  = [1 chirp_on];
chirp_off = [chirp_off size(chirp_offsets_temp,2)];

chirp_onsets  = chirp_onsets_temp(chirp_on);
chirp_offsets = chirp_offsets_temp(chirp_off);


%debug figure
% figure;
% plot(StimOn); ylim([-.5 1.5]);
% vline(stim_onsets(burst_on),'r')
% vline(stim_offsets(burst_off),'k')


ids.StimBounds{1}          = [intersect(chirp_onsets, find(Eye1On_only))' intersect(chirp_offsets, find(Eye1On_only))'];
ids.StimBounds{2}          = [intersect(chirp_onsets, find(Eye2On_only))' intersect(chirp_offsets, find(Eye2On_only))'];
ids.StimBounds{3}          = [intersect(chirp_onsets, find(bino_clean))' intersect(chirp_offsets, find(bino_clean))'];
ids.StimBoundsComment   = 'eye1, eye2, bino: [Chirp onset Chirp offset]';


% extract a few further stimparameters from aux_data
ids.trial_duration  = median(diff(frame_times_level(chirp_onsets))) / aux_samplingrate;
ids.stim_duration   = median(diff(frame_times_level([chirp_onsets' chirp_offsets']),[],2)) / aux_samplingrate;
ids.ITI_duration    = ids.trial_duration - ids.stim_duration;
ids.frame_times     = frame_times;
ids.level           = level;
ids.frame_times_level           = frame_times_level;
ids.aux_samplingrate = aux_samplingrate;