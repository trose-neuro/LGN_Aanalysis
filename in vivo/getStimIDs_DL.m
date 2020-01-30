function [ids varargout] =  getStimIDs_DL(aux_data, varargin)

% Not to self: THIS IS SUCH A HORRIBLE CODE!!! REWRITE ASAP
% PC = getenv('computername');

frame_cutoff = 0;
special = 0;
bscope2 = 1;
stimshutter = 1;

if ~isempty(varargin)
    numIndex = find(cellfun('isclass', varargin(1:end-1), 'char'));
    for ind = 1:length(numIndex)
        switch lower(varargin{numIndex(ind)})
            case 'eyes'
                eye = varargin{numIndex(ind) + 1}; %you never know when we are going to image this:http://goo.gl/y9pOH
            case 'level'
                level = varargin{numIndex(ind) + 1};
            case 'stimfile'
                stimfile = varargin{numIndex(ind) + 1};
            case 'aux_ch_frame_times'
                aux_ch_frame_times = varargin{numIndex(ind) + 1};
            case 'aux_ch_eye_times'
                aux_ch_frame_times = varargin{numIndex(ind) + 1};
            case 'aux_ch_stim'
                aux_ch_stim = varargin{numIndex(ind) + 1};
            case 'frame_cutoff'
                frame_cutoff = varargin{numIndex(ind) + 1};
            case 'special'
                special = varargin{numIndex(ind) + 1};
            case 'bscope'
                bscope2 = varargin{numIndex(ind) + 1};
            case 'stimshutter'
                stimshutter = varargin{numIndex(ind) + 1};
        end
    end
end



if bscope2
    aux_ch_frame_times = 4;
    aux_ch_eye_times = 18;
    aux_ch_eye_times2 = 17; %second eye rec will be inverted simply
%     aux_ch_eye_times2 =aux_ch_eye_times2 *-1 + max(aux_ch_eye_times2);
    aux_ch_stim = 8;
    aux_ch_stimshutter = 3;
    ids = struct;
    ids_new = struct;
else
    aux_ch_frame_times = 3;
    aux_ch_eye_times = 9;
    aux_ch_eye_times2 = 10;
    aux_ch_stim = 8;
    ids = struct;
    ids_new = struct;
end


if bscope2
     disp('USING NEW AUXDATA CHANNEL CONVENTIONS FOR BSCOPE2!');
end



% aux_data = load_lvd(['D:\data\mouse\sparse\' animal '\Data\' expdate filesep auxID]);
% stiminfo =  readini(stimfile);

if special
    disp('special case for messed up aux')
    aux_data(8,end-8000:end-7000) = 0;
    aux_data(8,end-7000:end-2000) = 5;
    aux_data(8,end-2000:end) = 0;
end

if bscope2
    aux_data(aux_ch_eye_times2,:) = aux_data(aux_ch_eye_times2,:)*-1+max(aux_data(aux_ch_eye_times2,:));
end

% stiminfo.seq =  repmat(stiminfo.seq,1,stiminfo.repetitions);
if max(stiminfo.seq) > 360
    disp('I assume 90 degree rotation (SCREEN PORTRAIT) and will subtract that from all orientations!')
    stiminfo.seq = stiminfo.seq - 90;
end

frame_times = get_frame_times(aux_data(aux_ch_frame_times,:));
%TR12 mod for multilevel ACQ
frame_times = frame_times(level:level:end); %TR12 for max stacks use the last frametime.

frame_times_a = frame_times;

%Be aware: depending on your experiment, the eye variable names might be
%wrong!

%first get the stimbins from the raw (not eye-specific) data.
if ~stimshutter
    stim_off = aux_data(aux_ch_stim, frame_times) - min(aux_data(aux_ch_stim, frame_times))< 0.1;
else
    stim_off = aux_data(aux_ch_stim, frame_times) - min(aux_data(aux_ch_stim, frame_times))< 0.1 & aux_data(aux_ch_stimshutter, frame_times) - min(aux_data(aux_ch_stimshutter, frame_times)) > 0.75;
end
stim_on_bin = aux_data(aux_ch_stim, frame_times) - min(aux_data(aux_ch_stim, frame_times)) > 0.75;

% get the stim off frame index
c_eye = sort([strfind(stim_on_bin,[1,0])]); %stim-binned frametimes
% c_eye3 = sort([strfind(stim_off+stim_on_bin,[0,1])]); %stim-binned frametimes

%account for random initial voltage in stimplotter
% if special == 1;
%     c_eye(1) = [];
% end

if eye == 2;
    %get the frametimes where one of the eyes is open
    eye_times =  aux_data(aux_ch_eye_times, frame_times) - min(aux_data(aux_ch_eye_times, frame_times)); %zero offset eye times in aux frametime samples
    eye2_times =  aux_data(aux_ch_eye_times2, frame_times) - min(aux_data(aux_ch_eye_times2, frame_times));
    %     eye_open{1} = eye_times < 0.1;
    eye_open{1} = eye_times > 0.75; % logical array for open eye 1 in aux frametimes
    eye_open{2} = eye2_times > 0.75;
    eye_frames{1} = find( eye_open{1} ==1);
    eye_frames{2} = find( eye_open{2} ==1);
    stimeye{1} = (eye_times(c_eye) > 0.75)'; % logical array denoting open eye 1 at start of new prestim period
    stimeye{2} = (eye2_times(c_eye) > 0.75)';
elseif eye == 3;
    %get the frametimes where one of the eyes or both eyes is open
    eye_times =  aux_data(aux_ch_eye_times, frame_times) - min(aux_data(aux_ch_eye_times, frame_times)); %eye  times in frametime
    eye2_times =  aux_data(aux_ch_eye_times2, frame_times) - min(aux_data(aux_ch_eye_times2, frame_times));
    eye3_times = eye_times + eye2_times;
    %     eye_open{1} = eye_times < 0.1;
    eye_open{1} = eye_times > 0.75 & eye3_times < 7;
    eye_open{2} = eye2_times > 0.75 & eye3_times < 7;
    eye_open{3} = eye3_times > 5.75;
    
    eye_frames{1} = find( eye_open{1} ==1);
    eye_frames{2} = find( eye_open{2} ==1);
    eye_frames{3} = find( eye_open{3} ==1);
    
    stimeye{1} = (eye_times(c_eye) > 0.75 & eye3_times(c_eye) < 7)';
    stimeye{2} = (eye2_times(c_eye) > 0.75 & eye3_times(c_eye) < 7)';
    stimeye{3} = (eye3_times(c_eye) > 5.75 )';
else
    eye_times = aux_data(aux_ch_eye_times, frame_times) - min(aux_data(aux_ch_eye_times, frame_times)); %eye  times in frametime;
    eye_open{1} = logical(ones(size(frame_times)));
    stimeye{1} = logical(ones(size(c_eye))');
end

% var_level = 0.1;% try to get the stim conditions from the analog signal

for eyes = 1:eye
    ids(eyes).eye_open = eye_open{eyes};
    
    frame_times = frame_times_a(eye_open{eyes}==1);
    
    ids(eyes).frame_times = frame_times;
    
    ids(eyes).stim_frames = aux_data(aux_ch_stim, frame_times) - min(aux_data(aux_ch_stim, frame_times));
    
    
    % stim_conditions = unique(fix(stim_frames * 10)) / 10;
    % stim_conditions = stim_conditions(stim_conditions > 0.75); % remove 0 and optional offset at the beginning
    % stim_on = uint8(zeros(1, length(stim_frames)));
    % for ind = 1:length(stim_conditions)
    %     stim_on(find((stim_frames > stim_conditions(ind) - var_level) & (stim_frames < stim_conditions(ind) + var_level))) = ind;
    % end
    %     ids(eyes).stim_off = ids(eyes).stim_frames < 0.1;
    
    if ~stimshutter
        ids(eyes).stim_off = ids(eyes).stim_frames < 0.1;
    else
        ids(eyes).stimshutterframes = aux_data(aux_ch_stimshutter, frame_times) - min(aux_data(aux_ch_stimshutter, frame_times));
        
        ids(eyes).stim_off = ids(eyes).stim_frames < 0.1 & ids(eyes).stimshutterframes > 0.75;
    end
    
    ids(eyes).stim_on_bin = ids(eyes).stim_frames > 0.75;
    c = sort([0 strfind(ids(eyes).stim_on_bin,[1,0]) length(ids(eyes).stim_on_bin)]);  %c -> Stimonsets including 0 and the last frame
    c = c(1:(stiminfo.orientations /eye * stiminfo.repetitions)+1); %truncate to stim presentations +1
    
    stimseq = (stiminfo.seq(stimeye{eyes}) / (360/length(unique(stiminfo.seq(stimeye{eyes}))))) + 1; %convert stim orientations to linear indices
    ids(eyes).stimseq_deg = stiminfo.seq(stimeye{eyes});
    ids(eyes).stim_on = [];
    for ii=1:length(c)-1,
        ids(eyes).stim_on = [ids(eyes).stim_on ids(eyes).stim_on_bin(c(ii)+1:c(ii+1))*stimseq(ii)]; %ids(eyes).stim_on is the linear stimulus code in frametime
    end
    ids(eyes).stim_on_deg = [];
    for ii=1:length(c)-1,
        temp = int16(ids(eyes).stim_on_bin(c(ii)+1:c(ii+1)));
        temp(temp == 0) = -1;
        temp(temp == 1) = ids(eyes).stimseq_deg(ii);
        ids(eyes).stim_on_deg = [ids(eyes).stim_on_deg temp]; %ids(eyes).stim_on_deg is the stimulus degree in frametime
    end
    
    ids(eyes).stim_on = [ids(eyes).stim_on zeros(1,length(ids(eyes).stim_on_bin) - length(ids(eyes).stim_on))]; %zeropad to length of stimonbin
    ids(eyes).stim_on_deg = [ids(eyes).stim_on_deg ones(1,length(ids(eyes).stim_on_bin) - length(ids(eyes).stim_on_deg)) * -1];
    stim_repetitions = sum(diff(ids(eyes).stim_on) == 1);
    
    if ~stimshutter
        bsline = [ 1 strfind(ids(eyes).stim_off,[0,1])];
        bsline = bsline(1:end-1);
        resp = strfind(ids(eyes).stim_on_bin,[0,1]);
        respend = strfind(ids(eyes).stim_on_bin,[1,0]);
    else
        bsline = strfind(ids(eyes).stim_off,[0,1]);
        bsline = bsline(1:2:end);
        resp = strfind(ids(eyes).stim_on_bin,[0,1]);
        respend = strfind(ids(eyes).stim_on_bin,[1,0]);
    end
    
    co = sort([bsline + 1 resp+1 respend length(ids(eyes).stim_on_bin)]); %co -> sorted combined StimON and StimOFFsets including 0 and the last frame
    
    %     co = co(1:(stiminfo.orientations /eye* stiminfo.repetitions * 2)+1);
    
    
    e = [];
    e(2:3:stiminfo.orientations /eye* stiminfo.repetitions * 3) = co(2:3:end); %stimulus start
    e(3:3:end+1) = co(3:3:end) - 1; %stimulus end
    e(1:3:end) = co(1:3:end-1); %baseline start
    %
    %     try
    %         e(1) = max(find((ids(eyes).stim_on_bin+ids(eyes).stim_off) == 0)) + 1;
    %     end
    
    ids(eyes).stim_boundaries = reshape(e, 3, length(e)/3)'; % stim baseline, start and end indices
    
    
    % let's reshuffle the stim to bring in order
    % this takes the actual stimorder and sorts the stimboundaries
    % according to increasing orientation angles in the rows and
    % repititions in the collums
    [~,f] = sort(reshape(stimseq, max(stimseq), stim_repetitions)',2); %#ok<UDIM>
    f = f + repmat(0:max(stimseq):max(stimseq)*(stim_repetitions - 1), max(stimseq), 1)';
    ids(eyes).stim_boundaries = ids(eyes).stim_boundaries(reshape(f',1,numel(f)),:); % reshuffle
    ids(eyes).stim_boundaries = reshape(ids(eyes).stim_boundaries, max(stimseq), stim_repetitions, 3);
    % stim_boundaries(stim category, repetitions, [stim_off_start stim_on_start stim_on_end])
    
    ids(eyes).stim_off_cond = zeros(size(ids(eyes).stim_off));
    ids(eyes).stim_on_bin_w_cutoff = logical(zeros(size(ids(eyes).stim_on)));
    
    for ind = 1:max(stimseq)
        
        for knd = 1:stim_repetitions
            ids(eyes).stim_off_cond(ids(eyes).stim_boundaries(ind,knd,1):ids(eyes).stim_boundaries(ind,knd,2)-1) = ind;
            ids(eyes).stim_on_bin_w_cutoff(ids(eyes).stim_boundaries(ind,knd,2) + frame_cutoff:ids(eyes).stim_boundaries(ind,knd,3)) = 1;
        end
        
    end
    
    
end

%% alternative second run (I'm lazy - I simply copy the code from above) to get not eye-divided timestamps for
% PSTHs

frame_times =frame_times_a;

for eyes = 1:eye;
    
    eye_open_new{1} = logical(ones(size(frame_times)));
    stimeye_new{1} = logical(ones(size(c_eye))');
    
    ids_new(eyes).eye_open_new = eye_open_new{1};
    
    frame_times = frame_times_a(eye_open_new{1}==1);
    
    ids_new(eyes).frame_times = frame_times;
    
    ids_new(eyes).stim_frames = aux_data(aux_ch_stim, frame_times) - min(aux_data(aux_ch_stim, frame_times));
    
    
    % stim_conditions = unique(fix(stim_frames * 10)) / 10;
    % stim_conditions = stim_conditions(stim_conditions > 0.75); % remove 0 and optional offset at the beginning
    % stim_on = uint8(zeros(1, length(stim_frames)));
    % for ind = 1:length(stim_conditions)
    %     stim_on(find((stim_frames > stim_conditions(ind) - var_level) & (stim_frames < stim_conditions(ind) + var_level))) = ind;
    % end
    %     ids_new(eyes).stim_off = ids_new(eyes).stim_frames < 0.1;
    
    if ~stimshutter
        ids_new(eyes).stim_off = ids_new(eyes).stim_frames < 0.1;
    else
        ids_new(eyes).stimshutterframes = aux_data(aux_ch_stimshutter, frame_times) - min(aux_data(aux_ch_stimshutter, frame_times));
        
        ids_new(eyes).stim_off = ids_new(eyes).stim_frames < 0.1 & ids_new(eyes).stimshutterframes > 0.75;
    end
    
    
    ids_new(eyes).stim_on_bin = ids_new(eyes).stim_frames > 0.75;
    
    c = sort([0 strfind(ids_new(eyes).stim_on_bin,[1,0]) length(ids_new(eyes).stim_on_bin)]);  %c -> Stimonsets including 0 and the last frame
    c = c(1:(stiminfo.orientations  * stiminfo.repetitions)+1); %truncate to stim presentations +1
    
    stimseq = (stiminfo.seq(stimeye_new{1}) / (360/length(unique(stiminfo.seq(stimeye_new{1})))) + 1); %convert stim orientations to linear indices
    ids_new(eyes).stimseq_deg = stiminfo.seq(stimeye_new{1});
    ids_new(eyes).stim_on = [];
    for ii=1:length(c)-1,
        ids_new(eyes).stim_on = [ids_new(eyes).stim_on ids_new(eyes).stim_on_bin(c(ii)+1:c(ii+1))*stimseq(ii)]; %ids_new(eyes).stim_on is the linear stimulus code in frametime
    end
    ids_new(eyes).stim_on_deg = [];
    for ii=1:length(c)-1,
        temp = int16(ids_new(eyes).stim_on_bin(c(ii)+1:c(ii+1)));
        temp(temp == 0) = -1;
        temp(temp == 1) = ids_new(eyes).stimseq_deg(ii);
        ids_new(eyes).stim_on_deg = [ids_new(eyes).stim_on_deg temp]; %ids_new(eyes).stim_on_deg is the stimulus degree in frametime
    end
    
    ids_new(eyes).stim_on = [ids_new(eyes).stim_on zeros(1,length(ids_new(eyes).stim_on_bin) - length(ids_new(eyes).stim_on))]; %zeropad to length of stimonbin
    ids_new(eyes).stim_on_deg = [ids_new(eyes).stim_on_deg ones(1,length(ids_new(eyes).stim_on_bin) - length(ids_new(eyes).stim_on_deg)) * -1];
    stim_repetitions = sum(diff(ids_new(eyes).stim_on) == 1)/eye;
    
    if ~stimshutter
        bsline = [ 1 strfind(ids_new(eyes).stim_off,[0,1])];
        bsline = bsline(1:end-1);
        resp = strfind(ids_new(eyes).stim_on_bin,[0,1]);
        respend = strfind(ids_new(eyes).stim_on_bin,[1,0]);
    else
        bsline = strfind(ids_new(eyes).stim_off,[0,1]);
        bsline = bsline(1:2:end);
        resp = strfind(ids_new(eyes).stim_on_bin,[0,1]);
        respend = strfind(ids_new(eyes).stim_on_bin,[1,0]);
    end
    
    co = sort([bsline + 1 resp+1 respend length(ids_new(eyes).stim_on_bin)]); %co -> sorted combined StimON and StimOFFsets including 0 and the last frame
    
    
    
%     co = sort([1 strfind(ids_new(eyes).stim_off,[1,0])+1 strfind(ids_new(eyes).stim_on_bin,[1,0])+1 length(ids_new(eyes).stim_on_bin)]); %co -> sorted combined StimON and StimOFFsets including 0 and the last frame
    
    %     co = co(1:(stiminfo.orientations * stiminfo.repetitions * 2)+1);
    
    onsets = co(2:3:end);
    offsets =  co(3:3:end) - 1;
    baseline = co(1:3:end-1);
    
    
    e=[];
    e(2:3:stiminfo.orientations / eye * stiminfo.repetitions * 3) = onsets(stimeye{eyes});
    e(3:3:end+1) = offsets(stimeye{eyes});
    e(1:3:end) = baseline(stimeye{eyes});
    %
    %     try
    %         e(1) = max(find((ids_new(eyes).stim_on_bin+ids_new(eyes).stim_off) == 0)) + 1;
    %     end
    
    ids_new(eyes).stim_boundaries = reshape(e, 3, length(e)/3)';
    
    
    % let's reshuffle the stim to bring in order
    % this takes the actual stimorder and sorts the stimboundaries
    % according to increasing orientation angles in the rows and
    % repititions in the collums
    [~,f] = sort(reshape(stimseq(stimeye{eyes}), max(stimseq(stimeye{eyes})), stim_repetitions)',2); %#ok<UDIM>
    f = f + repmat(0:max(stimseq(stimeye{eyes})):max(stimseq(stimeye{eyes}))*(stim_repetitions - 1), max(stimseq(stimeye{eyes})), 1)';
    ids_new(eyes).stim_boundaries = ids_new(eyes).stim_boundaries(reshape(f',1,numel(f)),:); % reshuffle
    ids_new(eyes).stim_boundaries = reshape(ids_new(eyes).stim_boundaries, max(stimseq(stimeye{eyes})), stim_repetitions, 3);
    % stim_boundaries(stim category, repetitions, [stim_off_start stim_on_start stim_on_end])
    
    ids_new(eyes).stim_off_cond = zeros(size(ids_new(eyes).stim_off));
    ids_new(eyes).stim_on_bin_w_cutoff = logical(zeros(size(ids_new(eyes).stim_on)));
    
    for ind = 1:max(stimseq(stimeye{eyes}))
        
        for knd = 1:stim_repetitions
            ids_new(eyes).stim_off_cond(ids_new(eyes).stim_boundaries(ind,knd,1):ids_new(eyes).stim_boundaries(ind,knd,2)-1) = ind;
            ids_new(eyes).stim_on_bin_w_cutoff(ids_new(eyes).stim_boundaries(ind,knd,2) + frame_cutoff:ids_new(eyes).stim_boundaries(ind,knd,3)) = 1;
        end
        
    end
    ids_new(eyes).stimseq_deg = stiminfo.seq(stimeye{eyes});
    ids_new(eyes).eye_open = ids_new(eyes).eye_open_new;
end

varargout{1} = ids_new;

