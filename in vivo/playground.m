roi_data = 'I:\David Laubender\Data\adata';
data  = 'I:\David Laubender\Data\imaging data';
savedir = 'C:\temp\figs'

mouse       = 'DL_191106_2';
% experiment  = 62335;
experiment  = 62336;


% aux channel settings
stim_channel = 8;
eye_channel  = 16;
bscope2 = 1;

% stim movie settings
im_size = 128; %pix
aspec   = 1.25;
dims = [im_size im_size*aspec];
getvids = 0;

adata_pathname      = fullfile(roi_data,mouse);
auxdata_pathname    = fullfile(data,mouse);

adata_file  = dir(fullfile(adata_pathname, ['**\*Adata*' num2str(experiment) '*.*']));
aux_file    = dir(fullfile(auxdata_pathname, ['**\*' num2str(experiment) '*.lvd']));


[eye2_ipsi eye1_contra stimvideo i1 i2 v1] = getvideos(num2str(experiment), aux_file.folder,3);



roi_data                           = load(fullfile(adata_file.folder, adata_file.name));
[auxdata aux_samplingrate]      = load_lvd(fullfile(aux_file.folder, aux_file.name));

level = str2num(cell2mat(regexp(roi_data.info.ImageDescription, '(?<=stackNumSlices = )\d+\.?\d*' , 'match' )));

ids = GET_StimIDs_Chirp( auxdata, aux_samplingrate, level )

if getvids
    frametimes = ids.frame_times;
    [vid_frame_times  sr]  = get_movie_frame_times(i1, auxdata,[], bscope2);
    [vid_frame_times2  sr] = get_movie_frame_times(i2, auxdata,[], bscope2);
    [svid_frame_times  sr] = get_movie_frame_times(v1, auxdata,[], bscope2);
    
    e1 = imresize(eye2_ipsi(:,:,vid_frame_times(1:level:end)),[dims(1) dims(2)]);
    e2 = imresize(eye1_contra(:,:,vid_frame_times2(1:level:end)),[dims(1) dims(2)]);
    v1 = imresize(stimvideo(:,:,svid_frame_times(1:level:end)),[dims(1) dims(2)]);
    stimstack = ([e1;e2;v1]);
    
    k =  view_tiff_stack_tr(stimstack);
end

figure(46363)

subplot(3,1,1), plot(smooth(roi_data.ROIs(140).activity))
subplot(3,1,2), plot(auxdata(8,ids.frame_times_level))
subplot(3,1,2), plot(auxdata(8,ids.frame_times_level))



%% PSTHing
eyes   = size(ids.StimBounds,2);
trials = size(ids.StimBounds{1},1);
rois   = size(roi_data.ROIs,2);

ratio = 0;
single_roi = 1;
plotdata = 1;
npfct = 0.7;
roi_data_cell{1} = roi_data;

for ROI = 1:rois
    PSTH_plot_Chirp(roi_data_cell, [], auxdata, ids, ROI, ratio, single_roi, plotdata, npfct, eyes, 0);
    saveas(gcf, [savedir 'psth_' num2str(ROI) '.png']);
end

% prestim  = floor(ids.ITI_duration / 2);
% poststim = floor(ids.ITI_duration / 2);
% stimdur  = ids.stim_duration;
% stimcutoff = 0.1;
%
% F0per = 1;
% detrend = 1;
% useF0mininterp = 1;
% df2nd           = 1;     % trace-by-trace secondary DF correction
% noneg           = 1;     % cut off at 0 ;




for roinum = 1:size(roi_data.ROIs,2)
    figure(roinum);
    for eyenum = 1:eyes
        for trialnum = 1:trials
            
        end
    end
end
