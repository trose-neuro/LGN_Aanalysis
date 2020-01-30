function        [DFpeak F0 F02 DFmean_trials SigF0 delta_max delta_mean delta_integ delta_meanF0 delta_sigmaF0 mean_plotdata plot_data_all] =  PSTH_plot_Puff(roi_data, ana_data, ids, ROI, ratio, single_roi, plotdata, npfct, eyes, basel, runexclude, varargin)

% disp(['extracting PSTH data of rec: ' ana_data{1}.info.FileID])
% plotdata = 0;
global roinum

suppress_disp = 1;

try
    if ROI == 1 && isempty(roinum);
        roinum = 1;
    end
end

nofit = 0;
try
    if nofit
        ana_data{1} = rmfield(ana_data{1},'Fit')
    end
end
recs = length(roi_data);

if isempty(npfct)
    npfct = 0.7;
end

if isempty(varargin)
    ROIimgs_ovl = []; %!!!ONLY RUN FOR VERY FIRST ANALYSIS OF ROIS!
    %     ROIimgs_green = [];
    ROIimgs_odi = [];
    
else
    ROIimgs_ovl = varargin{1};
    %     ROIimgs_green = varargin{2};
    ROIimgs_odi = varargin{2};
    
end

npfct_extract = 0;

pubfig = 1;
level = 4;

% plotdata = 1;

np_subtract = 1;

if ~np_subtract || npfct==0 && ~suppress_disp
    disp('No Neuropil Subtraction!')
end

poststim = ids{1}.ITI_duration / 4;
prestim = ids{1}.ITI_duration / 4;
stimdur = ids{1}.stim_duration;
stimcutoff = 0.1;
F0per = 1;

episodic_F0 = 0;
detrend = 1;
useF0mode = 0;
useF0mininterp = 1;

df2nd           = 1;     % trace-by-trace secondary DF correction
noneg           = 1;     % cut off at 0 ;

ylimover = 1;
ylo = [-200 1000];

% NP parameters
medianadd = 1;

endper = 1;
endsub = 0;
startper = .5;
timewidth = 0.5; %timewidth of filter in seconds

% bleedthrough corr
bleedcorr = 0;
if ~suppress_disp && ~bleedcorr
    disp('NOT correcting for g->R bleedthrough');
end
green2red_bleed = 0.05;


% if runexclude
dist_crit = 1;  %distance covered (a.u.) to count as running trial
speed_crit = 0.23; %2cm/s
mintrials = 4; % number of trials to be drawn after runexclusion (also true for anesth!)
% run_speed_thresh
% end

try
    if isfield(ana_data{1}, 'Fit')
        furtherfig = 4;
    else
        furtherfig = 2;
    end
end

expave = 0;
movave = 0;

downsamp = 1;

miny = -50;

minx = -1;

scalebarx = 10;
scalebary = 200;


if plotdata
    try
        close(23553);
        hdl = figure(23553);
    catch
        hdl = figure(23553);
    end
    
    %     figure
    pos = get(0, 'ScreenSize');
    %     pos(4) = (pos(4)-111)/2;
    %     pos(2) = 45+pos(4);
    %     set(gcf, 'Position' ,pos,  'Color', 'n');
    set(gcf, 'Position' ,pos);
    
    hold off;
end


k = 1;
recs =1;


for rec = 1:recs;
    
    SamplingFreq1 = regexp(ana_data{rec}.info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
    SamplingFreq1 = str2num(SamplingFreq1{1}) / level;
    
    %     smootho = round(timewidth * SamplingFreq1);
    %     cutoff = (0.44294/(sqrt(smootho^2-1)))*SamplingFreq1;
    smootho = 1;
    % disp(['moving average 3db cutoff: ' num2str(cutoff) ' Hz']);
    
    response = round(prestim*SamplingFreq1):round((prestim+stimdur-stimcutoff)*SamplingFreq1);
    baseline = round((prestim-F0per)*SamplingFreq1):round((prestim)*SamplingFreq1);
    endline = round((endper)*SamplingFreq1);
    startline = round((startper)*SamplingFreq1);
    
    
    if df2nd
        if roi_data{rec}.info.darkframes <100
            Dframes  = [roi_data{rec}.info.darkframes([1])+2:roi_data{rec}.info.darkframes(end)-2];
        else
            Dframes = [ceil(roi_data{rec}.info.darkframes(1)/4)+3: floor(roi_data{rec}.info.darkframes(end)/4)-3];
        end
        DFend    = nanmean(roi_data{rec}.ROIs(ROI(rec)).activity(Dframes));
        DFendr    = nanmean(roi_data{rec}.ROIs(ROI(rec)).activity_r(Dframes));
        DFnp      = nanmean(roi_data{rec}.np(ROI(rec)).activity(Dframes));
        DFnpr    = nanmean(roi_data{rec}.np(ROI(rec)).activity_r(Dframes));
        
        roi_data{rec}.ROIs(ROI(rec)).activity = roi_data{rec}.ROIs(ROI(rec)).activity - DFend;
        roi_data{rec}.ROIs(ROI(rec)).activity_r = roi_data{rec}.ROIs(ROI(rec)).activity_r - DFendr;
        roi_data{rec}.np(ROI(rec)).activity = roi_data{rec}.np(ROI(rec)).activity - DFnp;
        roi_data{rec}.np(ROI(rec)).activity_r = roi_data{rec}.np(ROI(rec)).activity_r - DFnpr;
    end
    
    if bleedcorr
        roi_data{rec}.ROIs(rec).activity_r = roi_data{rec}.ROIs(rec).activity_r - green2red_bleed * roi_data{rec}.ROIs(rec).activity;
        roi_data{rec}.np(rec).activity_r = roi_data{rec}.np(rec).activity_r - green2red_bleed * roi_data{rec}.np(rec).activity;
    end
    
    
    for ey = 1:eyes;
        
        ROIs = length(roi_data{rec}.ROIs);
        reps = size(ids{1}.StimBounds,1);
        
        plot_data = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));
        
        raw_data = roi_data{rec}.ROIs(ROI(rec)).activity;
        
        
        try
            np_data = roi_data{rec}.np(ROI(rec)).activity;
        catch
            %             disp('NO NEUROPIL! Using flatline!');
            roi_data{rec}.np(ROI(rec)).activity = zeros(size(full_raw_data));
            np_data = roi_data{rec}.np(ROI(rec)).activity;
        end
        
        if ~np_subtract
            roi_data{rec}.np(ROI(rec)).activity = zeros(size(raw_data));
            np_data = roi_data{rec}.np(ROI(rec)).activity;
        end
        
        full_np_data = roi_data{rec}.np(ROI(rec)).activity;
        
        if ratio
            plot_data_r = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));
        end
        
        if ratio
            raw_data_r = roi_data{rec}.ROIs(ROI(rec)).activity_r;
            np_data_r = roi_data{rec}.np(ROI(rec)).activity_r;
        end
        
        % neuropil correction (including median of np_trace
        % addition!)
        
        if npfct_extract
            %estimating the linear neuropil contamination factor similar to
            %the dendrite contribution in spine imaging data using
            %robustfit on the NP vs. cell samplewise plot.
            
            [int] = robustfit([np_data(200:end)], [raw_data(200:end)], 'andrews', [], 'off' );
            npfct = int;
            inter = 0;
            %
            %                         figure(21414)
            %                         plot(np_data(100:end), raw_data(100:end), '.k'); hold on;
            %                         plot(np_data(100:end), inter + npfct * np_data(100:end),'r','LineWidth',2);
            %                         xlabel('NP signal raw (F)')
            %                         ylabel('Cell signal raw (F)')
            disp(['slope ' num2str(npfct) 'intercept ' num2str(inter)]);
        else
            inter = 0;
        end
        
        if medianadd
            raw_data = ( raw_data - npfct * np_data) + npfct * nanmedian(np_data);
        else
            raw_data = ( raw_data - npfct * np_data);
        end
        
        if ratio
            if medianadd
                raw_data_r = (raw_data_r - npfct * np_data_r) + npfct * nanmedian(np_data_r);
            else
                raw_data_r = (raw_data_r - npfct * np_data_r);
            end
            % ratio
            raw_data = raw_data ./ raw_data_r;
        end
        
        if expave
            raw_data2 = tsmovavg(raw_data, 'e', smootho);
            raw_data = raw_data2;
            full_raw_data2 = tsmovavg(full_raw_data, 'e', smootho);
            full_raw_data = full_raw_data2;
        end
        
        if movave
            raw_data2 = smooth(raw_data,  smootho);
            raw_data = raw_data2;
            full_raw_data2 = smooth(full_raw_data,  smootho);
            full_raw_data = full_raw_data2;
        end
        %
        %         for ori = 1:oris;
        %             for rep = 1:reps
        %                 F0_new(ori,rep) = nanmedian(raw_data(ids{rec}(ey).stim_boundaries(ori,rep,2)-1-round(F0per*SamplingFreq1):ids{rec}(ey).stim_boundaries(ori,rep,2)-1)); %TR14 -1! before the first stimframe was included - that's of course wrong!
        %                 F0 = nanmedian(F0_new(:));
        %             end
        %         end
        %
        %
        
        if ~episodic_F0
            for rep = 1:reps
                baseperiodstim = ids{rec}(ey).StimBounds(rep,1)-1-round(F0per*SamplingFreq1):ids{rec}(ey).StimBounds(rep,1)-1;
                if min(baseperiodstim) <1
                    F0_new(rep) = NaN;
                else
                    
                    
                    F0_new(rep) = nanmedian(raw_data(baseperiodstim)); %TR14 -1! before the first stimframe was included - that's of course wrong!
                    F0 = nanmedian(F0_new(:));
                end
                if useF0mode
                    disp('Using mode of F0 distributions as global F0');
                    F0 = mode(F0_new(:));
                end
            end
        end
        
        
        if noneg;
            raw_data(raw_data <=0) = 0.01;
        end
        
        
        for rep = 1:reps;
            
            plotarray_disp =  ids{rec}(ey).StimBounds(rep,1)-1-round(prestim*SamplingFreq1):ids{rec}(ey).StimBounds(rep,2)+floor(SamplingFreq1*poststim);
            if min(plotarray_disp) <1
                plot_data(rep,:) = NaN;
                continue
            end
            % %
            %                         if ori == 6 && sf == 1 && tf == 1;
            %                             disp('hols');
            %                         end
            
            DF = raw_data(plotarray_disp)-F0;
            
            plot_data(rep,1:length(DF)) =(DF./F0) * 100;
            
            timebase = (1:length(plot_data)) / SamplingFreq1;
            
            %% consolidate data output fdrom single trace
            DFpeak(rep)         = nanmax(plot_data(rep,response));
            DFmean_trials(rep)  = nanmean(plot_data(rep,response)); % for ANOVA responder crit
            F02(rep)            = nanmean(plot_data(rep,baseline));
            SigF0(rep)          = nanstd(plot_data(rep,baseline));
            Fend(rep)           = nanmean(plot_data(rep,end-endline:end));
            Fstart(rep)         = nanmean(plot_data(rep,1:startline));
            
            
            if plotdata
                hdl(1) = subplot(1,1,k);
                
                ptr = plot(timebase(1:downsamp:end),plot_data(rep,1:downsamp:end), 'k');
                %                             try
                %                                 if  ids{rec}(ey).running_maxspeed(ori,rep) > speed_crit
                %                                     set(ptr, 'Color', 'b');
                %                                 end
                %                             end
                hold on;
            end
        end
        
        
        mean_plotdata{1} = nanmean(plot_data,1);
        plot_data_all{1} = plot_data;
        sum_plotdata{1}  = sum(plot_data,1);
        
        plot_data = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));  % cleaing up the incrementing plot_data array - otherwise old stims get carried over to the mean!
        
        timebase = (1:length(mean_plotdata{1})) / SamplingFreq1;
        
        if plotdata
            hdl(2) = subplot(1,1,k);
            plot(timebase(1:downsamp:end), mean_plotdata{1}(1:downsamp:end), 'r', 'LineWidth', 4);
            k = k+1;
            
            xlim([minx timebase(end)]);
            
            if~single_roi
                %                             lazy = [miny nanmax(displaypeaks(:))+1*nanmax(displaypeaks(:))];
                %
                %                             ylim(lazy);
                %                             if ylimover
                %                                 ylim(ylo)
                %                             end
                %
                %                             if ~pubfig;
                %                                 hline(displaypeaks(rec,ori,ey));
                %                             end
            else
                %                     hline(sum(sum_plotdata)/recs/stimdur);
                lazy = [0 1000];
                
                ylim(lazy);
                if ylimover
                    ylim(ylo)
                end
            end
            
            
            set(gca, 'visible', 'off');
            set(gca, 'Color', 'none');
            
            
            try
                if ana_data{rec}.peaks(ROI(rec)).responder_contra_n || ana_data{rec}.peaks(ROI(rec)).responder_ipsi_n
                    ylabel(['exp' ana_data{rec}.info.FileID ' * ']);
                else
                    ylabel(['exp' ana_data{rec}.info.FileID]);
                end
            end
            set(get(gca,'ylabel'),'Visible','on', 'FontSize', 14, 'Interpreter', 'none')
            
            %                         axis tight;
            
            yl = ylim;
            %                         yl = [miny miny+25];
            x = [prestim prestim+stimdur];
            y = [yl(2) yl(2)];
            y2 = [yl(1) yl(1)];
            
            patch([x fliplr(x)],[y fliplr(y2)], 'g' ,'facealpha',0.1,...
                'edgecolor','none', 'Tag', 'Stim');
            
            
            %make scalebars
            
            xs = [minx minx+1];
            xst = [minx minx+scalebarx];
            
            ys1 = [miny miny];
            ys2 = ys1+scalebary;
            
            yst1 = [miny miny];
            yst2  = [-4 -4];
            ptsb1 = patch([xs fliplr(xs)],[ys1 fliplr(ys2)], 'k',...
                'edgecolor','none', 'Tag', 'scalebary');
            
            ptsb2 = patch([xst fliplr(xst)],[yst1 fliplr(yst2)], 'k',...
                'edgecolor','none', 'Tag', 'scalebarx');
        end
        %% consolidate data output fdrom average trace
        
        if movave % if the full trace is smoothed, there is no reason to smooth the average trace
            smootho2 = 1;
        else
            smootho2 = smootho;
        end
        
        delta_max(1) = nanmax(smooth(mean_plotdata{1}(response),smootho2));
        delta_mean(1) = nanmean(smooth(mean_plotdata{1}(response),smootho2));
        delta_integ(1) = trapz(smooth(mean_plotdata{1}(response),smootho2));
        delta_meanF0(1) = nanmean(mean_plotdata{1}(baseline));
        delta_sigmaF0(1) = nanstd(mean_plotdata{1}(baseline));
        delta_meanrate(1) = sum(sum_plotdata{1}(response)) / reps / stimdur;
        
    end
end


if runexclude
    DFpeak =  DFpeak_re;
    F02 =  F02_re;
    SigF0 =  SigF0_re;
    Fend =  Fend_re;
    Fstart =  Fstart_re;
end

if plotdata
    if single_roi
        %                 set(hdl(:), 'ylim', [miny nanmax(delta_meanrate(:)) + 0.333*nanmax(delta_meanrate(:))]);
        set(hdl(:), 'ylim', [miny nanmax(delta_max(:)) + 0.333*nanmax(delta_max(:))]);
        if ylimover
            set(hdl(:), 'ylim', ylo);
        end
        
        tyl = ylim;
        if tyl(2)<miny+scalebary;
            set(hdl(:), 'ylim', [miny miny + scalebary + miny]);
            if ylimover
                set(hdl(:), 'ylim', ylo);
            end
        end
        tyl = ylim;
        txl = xlim;
        
        yl = ylim;
        y = [yl(2) yl(2)];
        y2 = [yl(1) yl(1)];
        p = findobj(gcf, 'Tag','Stim');
        set(p(:), 'ydata', [y fliplr(y2)]);
        
        p2 = findobj(gcf, 'Tag','scalebarx');
        p3 = findobj(gcf, 'Tag','scalebary');
        set(p2(:), 'ydata', [[miny miny] [miny+tyl(2)* 0.01 miny+tyl(2)* 0.01]]);
        set(p3(:), 'xdata', [[minx minx+txl(2)* 0.01 ] [minx+txl(2)* 0.01 minx]]);
        
        t = findobj(gcf, 'Type','text');
        poss = get(t(:), 'Position');
        for u = 1:length(poss);
            set(t(u), 'position', [poss{u}(1)+0.5  tyl(2) - 0.1 * tyl(2) 0],'FontSize',35)
        end
        hline(0,'--k');
    end
end


if plotdata
    tightfig;
    movegui('center');
end

roinum = roinum + 1;
% fillPage(gcf, 'papersize', 'A4'); orient landscape;
