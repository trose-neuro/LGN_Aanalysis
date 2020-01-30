function [F0_min_new F0_new F0_mini_interp F0_interp] = getF0frombaseline_chirp(raw_data, ids_temp, F0per, SamplingFreq1);
eyes = 2; %this for now only works for alternating monocular stimulation

try
    ids = ids{1};
end

eyes =size(ids_temp.StimBounds,2);
% simply reshape ids of chirp

ids(1).stim_boundaries(:,:,3) = ids_temp.StimBounds{1}(:,2);
ids(1).stim_boundaries(:,:,2) = ids_temp.StimBounds{1}(:,1);

ids(2).stim_boundaries(:,:,3) = ids_temp.StimBounds{2}(:,2);
ids(2).stim_boundaries(:,:,2) = ids_temp.StimBounds{2}(:,1);

if eyes == 3
    ids(3).stim_boundaries(:,:,3) = ids_temp.StimBounds{3}(:,2);
    ids(3).stim_boundaries(:,:,2) = ids_temp.StimBounds{3}(:,1);
end

%% concatenation of all presstimulus period indices (both eyes)
e1_base = [ids(1).stim_boundaries(:,:,2)-1-round(F0per*SamplingFreq1)];%:ids(1).stim_boundaries(:,:,2)-1;
e1_base = sort([e1_base(:)]);
e1_baseend = ids(1).stim_boundaries(:,:,2)-1;
e1_baseend = sort([e1_baseend(:)]);

e2_base = [ids(2).stim_boundaries(:,:,2)-1-round(F0per*SamplingFreq1)];%:ids(1).stim_boundaries(:,:,2)-1;
e2_base = sort([e2_base(:)]);
e2_baseend = ids(2).stim_boundaries(:,:,2)-1;
e2_baseend = sort([e2_baseend(:)]);

if eyes == 3;
    e3_base = [ids(3).stim_boundaries(:,:,2)-1-round(F0per*SamplingFreq1)];%:ids(1).stim_boundaries(:,:,2)-1;
    e3_base = sort([e3_base(:)]);
    e3_baseend = ids(3).stim_boundaries(:,:,2)-1;
    e3_baseend = sort([e3_baseend(:)]);
    
    base_idx(:,1) = sort([e1_base ; e2_base; e3_base]);
    base_idx(:,2) = sort([e1_baseend ; e2_baseend;e3_baseend ]);
else
    
    base_idx(:,1) = sort([e1_base ; e2_base]);
    base_idx(:,2) = sort([e1_baseend ; e2_baseend]);
    
end

%% concatenation of stimulus period indices (both eyes)
e1_stim = [ids(1).stim_boundaries(:,:,2)];%:ids(1).stim_boundaries(:,:,2)-1;
e1_stim = sort([e1_stim(:)]);
e1_stimend = ids(1).stim_boundaries(:,:,3);
e1_stimend = sort([e1_stimend(:)]);

e2_stim = [ids(2).stim_boundaries(:,:,2)];%:ids(1).stim_boundaries(:,:,2)-1;
e2_stim = sort([e2_stim(:)]);
e2_stimend = ids(2).stim_boundaries(:,:,3);
e2_stimend = sort([e2_stimend(:)]);


if eyes == 3;
    e3_stim = [ids(3).stim_boundaries(:,:,2)];%:ids(1).stim_boundaries(:,:,2)-1;
    e3_stim = sort([e3_stim(:)]);
    e3_stimend = ids(3).stim_boundaries(:,:,3);
    e3_stimend = sort([e3_stimend(:)]);
    
    stim_idx(:,1) = sort([e1_stim ; e2_stim; e3_stim]);
    stim_idx(:,2) = sort([e1_stimend ; e2_stimend; e3_stimend]);
else
    stim_idx(:,1) = sort([e1_stim ; e2_stim]);
    stim_idx(:,2) = sort([e1_stimend ; e2_stimend]);
    
end




%% extraction and subtraction of baseline ratios
raw_data_new = NaN(size(raw_data));

reclength = length(raw_data);
timebase = linspace(0, reclength/SamplingFreq1, reclength);

for i = 1:length(base_idx)
    raw_data_new(base_idx(i,1):base_idx(i,2)) = raw_data(base_idx(i,1):base_idx(i,2));
    
    F0(i) = nanmedian(raw_data(base_idx(i,1):base_idx(i,2)));
    
    x(i) = mean(base_idx(i,1):base_idx(i,2));
end

%% find minima and interpolate F0 trace.
[~,~,~, minidx] = extrema(F0);
minidx = unique([1 minidx length(F0)]); %just take the first and last F0 value to prevent stupid interpolation;

mini = F0(sort(minidx))';


F0_mini_interp  = interp1(x(sort(minidx)), smooth(mini,10,'rloess'), 1:reclength ,'linear', 'extrap');
%F0_interp       = interp1(x, smooth(F0,10,'rloess'), 1:reclength ,'linear', median(F0));
F0_interp = 1;
%% reextract F0 values from interpolated F0 traces.

for ey = 1:eyes;
    
    reps = size(ids(ey).stim_boundaries,1);
    
    for rep = 1:reps
        %F0_new(ori,rep,ey,ey) = nanmean(F0_interp(ids(ey).stim_boundaries(ori,rep,2)-1-round(F0per*SamplingFreq1):ids(ey).stim_boundaries(ori,rep,2)-1));
        F0_min_new(rep,ey) = nanmean(F0_mini_interp(ids(ey).stim_boundaries(rep,2)-1-round(F0per*SamplingFreq1):ids(ey).stim_boundaries(rep,2)-1));
    end
end
% figure(2352)
% plot(raw_data);hold all
% plot(F0_mini_interp); hold off
F0_new = 1;





