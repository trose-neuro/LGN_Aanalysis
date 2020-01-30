function ids = GET_StimIDs_SF_TF( AuxData, stimarray, level )


%  channels
Frames = AuxData(3,:);
Stims = AuxData(8,:);

frame_times = get_frame_times(Frames);
frame_times_level = frame_times(1:level:end);


StimOn = Stims(frame_times_level)>0.8;
stim_onsets = find(diff(StimOn==1)>0);
stim_offsets = find(diff(StimOn==1)<0);

if length(stim_offsets) < length(stim_onsets);
    %debug figure
    figure;
    plot(StimOn);
    ylim([-.5 1.5]);
    vline(stim_offsets(end), 'b');
    stim_offsets(length(stim_onsets)) = stim_offsets(end) + median(diff([stim_onsets stim_offsets]));
    vline(stim_offsets(end), 'r');
end

% Get full stored stimulus sequence

% Direction

stim_dirs = sort(unique(stimarray.save_angles(stimarray.repruns,:)));
stim_dirs_rand = stimarray.save_angles(stimarray.repruns,:);

if  size(stimarray.save_angles(stimarray.repruns,:),2) ~= size(stim_onsets,2);
    disp('Extracted stimuli number does ot fit stimulus file')
end

% spatial frequency
stim_SFs = sort(unique(stimarray.save_sf(stimarray.repruns,:)));
stim_SFs_rand = stimarray.save_sf(stimarray.repruns,:);

% temporal frequency
stim_TFs = sort(unique(stimarray.save_tf(stimarray.repruns,:)));
stim_TFs_rand = stimarray.save_tf(stimarray.repruns,:);


% Stimulus presentations

StimBounds = zeros(size(stim_dirs,2),2,size(stim_SFs,2),size(stim_TFs,2));

for Nr = 1:length(stim_onsets)
    Curr_Stim_dir = stim_dirs_rand(Nr);
    Curr_Stim_sf = stim_SFs_rand(Nr);
    Curr_Stim_tf = stim_TFs_rand(Nr);
    
    Stim_Number_dir = find(abs(stim_dirs - Curr_Stim_dir) == min(abs(stim_dirs - Curr_Stim_dir))); % index into dir stim
    Stim_Number_sf = find(abs(stim_SFs - Curr_Stim_sf) == min(abs(stim_SFs - Curr_Stim_sf))); % index into sf stim
    Stim_Number_tf = find(abs(stim_TFs - Curr_Stim_tf) == min(abs(stim_TFs - Curr_Stim_tf))); % index into tf stim
    
    trial = size(find(StimBounds(Stim_Number_dir,1,Stim_Number_sf,Stim_Number_tf,:)),1);
    if ~trial;
        StimBounds(Stim_Number_dir,1,Stim_Number_sf,Stim_Number_tf, 1) = stim_onsets(Nr);
        StimBounds(Stim_Number_dir,2,Stim_Number_sf,Stim_Number_tf, 1) = stim_offsets(Nr);
    else
        StimBounds(Stim_Number_dir,1,Stim_Number_sf,Stim_Number_tf, trial + 1 ) = stim_onsets(Nr);
        StimBounds(Stim_Number_dir,2,Stim_Number_sf,Stim_Number_tf, trial + 1 ) = stim_offsets(Nr);
    end
end

ids.StimBounds = StimBounds;
ids.StimBoundsComment = 'StimBounds Array Dimensions: DirectionOnsetFrame; DirectionOffsetFrame; SpatialFreq; TemporalFreq; Trial';
ids.stim_BLOCK = stimarray.repruns;
ids.stim_duration = stimarray.movieDurationSecs;
ids.stim_ITI = stimarray.inter_stim_interval;
ids.stim_dirs = stim_dirs;
ids.stim_SFs = stim_SFs;
ids.stim_TFs = stim_TFs;
ids.stim_dirs_seq = stim_dirs_rand;
ids.stim_SFs_seq = stim_SFs_rand;
ids.stim_TFs_seq = stim_TFs_rand;
ids.full_stimarray = stimarray;


%
% for s =  1:length(StimOnsets)
%     ids(eyes).StimOnsetFrames{s} = StimOnsetFrames{s}(StimSettings.Eyematrix==eyes-1);
%     ids(eyes).StimOffsetFrames{s} = StimOffsetFrames{s}(StimSettings.Eyematrix==eyes-1);
% end
%
% ids(eyes).eye_open = eye_open{eyes};
% frame_times = FrameOnsetsr(eye_open{eyes}==1);
% ids(eyes).frame_times = frame_times;
% ids(eyes).stim_frames = Stims(frame_times) - min(Stims(frame_times));
