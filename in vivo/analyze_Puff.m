function analyze_Puff;
%% Analyze Simon's airpuff data


%% load the single session .mat files and consolidate them
close all; clear all
adata_dir = 'I:\Simon Weiler\AnalyzedData'

cd(adata_dir);
s =rdir('**\*PUFF*.mat');

for cons = 1:length(s);
    concat_block(cons) = load(s(cons).name, '-mat');
    nROIs(cons) = length(concat_block(cons).peaks);
    nreps(cons) = size(concat_block(cons).ids.StimBounds,1);
end

%% extract a few stim parameters
poststim = concat_block(cons).ids.ITI_duration / 4;
prestim = concat_block(cons).ids.ITI_duration / 4 ;
stimdur = concat_block(cons).ids.stim_duration;
stimcutoff = 0;
endper = 1;
F0per = 1;
startper = .5;
SamplingFreq1 = concat_block(cons).info.SamplingFreq1;

plotlength = ceil((prestim+stimdur+poststim+stimcutoff+endper)*SamplingFreq1);
plot_data = NaN(1,plotlength);


%% consolidate everything in one structure, skipping excluded traces; calculating significance
k =1; p = 1; o = 1;
mean_concat_traces = [];
all_concat_traces = [];
mean_sig_concat_traces = [];
% whichrois = 19;

for recons = 1:length(concat_block);
    whichrois = 1:nROIs(recons);
    for rois = whichrois;
        if ~concat_block(recons).peaks(rois).excluded
            for reps =1:nreps(recons);
                all_concat_traces(k, :) = concat_block(recons).peaks(rois).plot_data_all{1}(reps,:);
                k = k + 1;
            end
            mean_concat_traces(p, :) = concat_block(recons).peaks(rois).mean_plotdata{1};
            mean_concat_meanPSTH(p, :) = concat_block(recons).peaks(rois).delta_mean;
            mean_concat_sigmaPSTHbaseline(p, :) = concat_block(recons).peaks(rois).delta_sigmaF0;            
            if concat_block(recons).peaks(rois).ANOVA < .05
                mean_sig_concat_traces(o, :) = concat_block(recons).peaks(rois).mean_plotdata{1};
                mean_sig_concat_meanPSTH(o, :) = concat_block(recons).peaks(rois).delta_mean;
                mean_concat_sigma_sig_PSTHbaseline(o, :) = concat_block(recons).peaks(rois).delta_sigmaF0;
                o = o + 1;
            end            
            p = p + 1;
        end
    end
end

%% baseline z-scoring
mean_concat_traces_z    = mean_concat_traces ./ repmat(mean_concat_sigmaPSTHbaseline, 1,size(mean_concat_traces,2));
mean_sig_concat_traces_z  = mean_sig_concat_traces ./ repmat(mean_concat_sigma_sig_PSTHbaseline,1,size(mean_sig_concat_traces,2));

%% timebase
time_base = (1:size(mean_sig_concat_traces,2)) / SamplingFreq1;

%% sorting for pixel display
[~, sortIDX] =  sort(mean_sig_concat_meanPSTH);
[~, sortIDX2] =  sort(mean_concat_meanPSTH);

%% make figure and extract the PSTH values

response = round(prestim*SamplingFreq1):round((prestim+stimdur-stimcutoff)*SamplingFreq1);
baseline = round((prestim-F0per)*SamplingFreq1):round((prestim)*SamplingFreq1);
endline = round((endper)*SamplingFreq1);
startline = round((startper)*SamplingFreq1);

% all traces
figure(3456356)
ll = shadedErrorBar(time_base ,mean_concat_traces_z,{@mean,@SEM},{'k','Marker','none', 'LineWidth', 2});
oo = vline([response(1)/SamplingFreq1 response(end)/SamplingFreq1]);
set(oo, 'LineWidth', 4)
xlabel('time')
ylabel('\DeltaF/F_0')
title('All traces')

view_tiff(mean_concat_traces(sortIDX2,:));
ppp = vline([response(1) response(end)]);
set(ppp, 'LineWidth', 4)
set(gcf, 'Name', 'all traces sorted')

view_tiff(mean_concat_traces_z(sortIDX2,:));
ppp = vline([response(1) response(end)]);
set(ppp, 'LineWidth', 4)
set(gcf, 'Name', 'all z-scored traces sorted')

% all significant traces
figure(34356356)
kk = shadedErrorBar(time_base ,mean_sig_concat_traces_z,{@mean,@SEM},{'k','Marker','none', 'LineWidth', 2})
pp = vline([response(1)/SamplingFreq1 response(end)/SamplingFreq1]);
set(pp, 'LineWidth', 4)
title('All significant traces')

view_tiff(mean_sig_concat_traces(sortIDX,:));
pppp = vline([response(1) response(end)]);
set(pppp,'LineWidth', 4)
set(gcf, 'Name', 'all significant traces sorted')

view_tiff(mean_sig_concat_traces_z(sortIDX,:));
pppp = vline([response(1) response(end)]);
set(pppp,'LineWidth', 4)
set(gcf, 'Name', 'all significant z-scored traces sorted')

disp('break');

%% functions
function serr = SEM(x);
serr = std(x)/sqrt(length(x));
