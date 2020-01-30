function        [DFpeak F0 F02 DFmean_trials SigF0 delta_max delta_mean delta_integ delta_meanF0 delta_sigmaF0 mean_plotdata plot_data_all] =  PSTH_plot_chronic(roi_data, ana_data, ids, ROI, ratio, single_roi, plotdata, npfct, eyes, basel, runexclude, varargin)

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

poststim = ids{1}.stim_ITI ;
prestim = ids{1}.stim_ITI ;
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


stims = size(ids{1}.stim_dirs,2);

if stims == 8
    x_sym_strings = 'VXJLNPRT';
    x_col = 'ygbrygbr';
elseif stims == 12
    x_sym_strings = 'VWYJKMNOQRSU';
    x_col = 'kkkkkkkkkkkk';
    ht = makeHueTable(12,1/12);
    ylix = find(ht(:,1) == 1 & ht(:,2) == 1 ); %quick hack: find yellow for vertical;
    x_col = makeHueTable(12,ylix/12);
elseif stims == 16
    x_sym_strings = 'VWXYJKLMNOPQRSTU';
    x_col = 'kkkkkkkkkkkkkkkk';
end

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

if ~single_roi
    for rec = 1:recs
        displaypeaks(rec,:,1) = [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_ipsi];
        displaypeaks(rec,:,2) = [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_contra];
        if eyes == 3;
            displaypeaks(rec,:,3) = [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_bino];
        end
    end
else
    recs =1;
end


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
        oris = size(ids{rec}(ey).stim_dirs,2);
        sfs = size(ids{rec}(ey).stim_SFs,2);
        tfs = size(ids{rec}(ey).stim_TFs,2);
        
        reps = size(ids{1}.StimBounds,5);
        
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
            for ori = 1:oris;
                for sf= 1:sfs;
                    for tf= 1:tfs;
                        % %             generate episodic prestim F0 data and take median of these for overall F0 value of detrended raw data
                        for rep = 1:reps
                            baseperiodstim = ids{rec}(ey).StimBounds(ori,1,sf,tf,rep)-1-round(F0per*SamplingFreq1):ids{rec}(ey).StimBounds(ori,1,sf,tf,rep)-1;
                            if min(baseperiodstim) <1
                                F0_new(ori,sf,tf,rep) = NaN;
                            else
                                
                                
                                F0_new(ori,sf,tf,rep) = nanmedian(raw_data(baseperiodstim)); %TR14 -1! before the first stimframe was included - that's of course wrong!
                                F0 = nanmedian(F0_new(:));
                            end
                            if useF0mode
                                disp('Using mode of F0 distributions as global F0');
                                F0 = mode(F0_new(:));
                            end
                        end
                    end
                    
                end
            end
        end
        
        if noneg;
            raw_data(raw_data <=0) = 0.01;
        end
        
        
        for  sf = 1:sfs;
            for tf = 1:tfs;
                for ori = 1:oris;
                    for rep = 1:reps;
                        
                        plotarray_disp =  ids{rec}(ey).StimBounds(ori,1,sf,tf,rep)-1-round(prestim*SamplingFreq1):ids{rec}(ey).StimBounds(ori,2,sf,tf,rep)+floor(SamplingFreq1*poststim);
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
                        DFpeak(ori,sf,tf,rep) = nanmax(plot_data(rep,response));
                        DFmean_trials(ori,sf,tf,rep) = nanmean(plot_data(rep,response)); % for ANOVA responder crit
                        F02(ori,sf,tf,rep) = nanmean(plot_data(rep,baseline));
                        SigF0(ori,sf,tf,rep) = nanstd(plot_data(rep,baseline));
                        Fend(ori,sf,tf,rep) = nanmean(plot_data(rep,end-endline:end));
                        Fstart(ori,sf,tf,rep) = nanmean(plot_data(rep,1:startline));
                        
                        
                        if plotdata
                            hdl(ori,sf,tf) = subplot(sfs, oris * tfs ,k);
                            
                            ptr = plot(timebase(1:downsamp:end),plot_data(rep,1:downsamp:end), 'k');
                            try
                                if  ids{rec}(ey).running_maxspeed(ori,rep) > speed_crit
                                    set(ptr, 'Color', 'b');
                                end
                            end
                            hold on;
                        end                       
                    end
                    
                    
                    mean_plotdata{ori,sf,tf} = nanmean(plot_data,1);
                    plot_data_all{ori,sf,tf} = plot_data;
                    sum_plotdata{ori,sf,tf} = sum(plot_data,1);

                    plot_data = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));  % cleaing up the incrementing plot_data array - otherwise old stims get carried over to the mean!
                    
                    timebase = (1:length(mean_plotdata{ori,sf,tf})) / SamplingFreq1;
                    
                    if plotdata
                        hdl(ori,sf,tf) = subplot(sfs, oris *tfs ,k);
                        plot(timebase(1:downsamp:end), mean_plotdata{ori,sf,tf}(1:downsamp:end), 'r', 'LineWidth', 1);
                        k = k+1;
                        
                        xlim([minx timebase(end)]);
                        
                        if~single_roi
                            lazy = [miny nanmax(displaypeaks(:))+1*nanmax(displaypeaks(:))];
                            
                            ylim(lazy);
                            if ylimover
                                ylim(ylo)
                            end
                            
                            if ~pubfig;
                                hline(displaypeaks(rec,ori,ey));
                            end
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
                        
                        if ori == 1 & ey ==1;
                            
                            try
                                if ana_data{rec}.peaks(ROI(rec)).responder_contra_n || ana_data{rec}.peaks(ROI(rec)).responder_ipsi_n
                                    ylabel(['exp' ana_data{rec}.info.FileID ' * ']);
                                else
                                    ylabel(['exp' ana_data{rec}.info.FileID]);
                                end
                            end
                            set(get(gca,'ylabel'),'Visible','on', 'FontSize', 14, 'Interpreter', 'none')
                        end
                        %                         axis tight;
                        
                        yl = ylim;
                        %                         yl = [miny miny+25];
                        x = [prestim prestim+stimdur];
                        y = [yl(2) yl(2)];
                        y2 = [yl(1) yl(1)];
                        
                        %                         if tf ==1 ;
                        %                             patch([x fliplr(x)],[y fliplr(y2)], [0.044 0.162 0.095] ,'facealpha',0.3,...
                        %                                 'edgecolor','none', 'Tag', 'Stim');
                        %                         elseif tf == 2 ;
                        %                             patch([x fliplr(x)],[y fliplr(y2)], [0.153 0.216 0.201] ,'facealpha',0.3,...
                        %                                 'edgecolor','none', 'Tag', 'Stim');
                        %                         elseif tf == 3;
                        %                             patch([x fliplr(x)],[y fliplr(y2)], [0.229 0.245 0.249] ,'facealpha',0.3,...
                        %                                 'edgecolor','none', 'Tag', 'Stim');
                        %                         end
                        %
                        if tf ==1 ;
                            patch([x fliplr(x)],[y fliplr(y2)], 'g' ,'facealpha',0.1,...
                                'edgecolor','none', 'Tag', 'Stim');
                        elseif tf == 2 ;
                            patch([x fliplr(x)],[y fliplr(y2)], 'b' ,'facealpha',0.1,...
                                'edgecolor','none', 'Tag', 'Stim');
                        elseif tf == 3;
                            patch([x fliplr(x)],[y fliplr(y2)], 'r','facealpha',0.1,...
                                'edgecolor','none', 'Tag', 'Stim');
                        end
                        
                        %make scalebars
                        if ori == 1 && ey ==1;
                            xs = [minx minx+1];
                            xst = [minx minx+scalebarx];
                            
                            ys1 = [miny miny];
                            ys2 = ys1+scalebary;
                            
                            yst1 = [miny miny];
                            yst2  = [-40 -40];
                            ptsb1 = patch([xs fliplr(xs)],[ys1 fliplr(ys2)], 'k',...
                                'edgecolor','none', 'Tag', 'scalebary');
                            
                            ptsb2 = patch([xst fliplr(xst)],[yst1 fliplr(yst2)], 'k',...
                                'edgecolor','none', 'Tag', 'scalebarx');
                            
                            %                                         disp(['scalebar: ' num2str(scalebarx) ' s - ' num2str(scalebary) '% DR/R0']);
                            
                        end
                        %                 ylabel(['exp' roi_data{rec}.fnames])
                    end
                    %% consolidate data output fdrom average trace
                    
                    if movave % if the full trace is smoothed, there is no reason to smooth the average trace
                        smootho2 = 1;
                    else
                        smootho2 = smootho;
                    end
                    
                    delta_max(ori,sf,tf) = nanmax(smooth(mean_plotdata{ori,sf,tf}(response),smootho2));
                    delta_mean(ori,sf,tf) = nanmean(smooth(mean_plotdata{ori,sf,tf}(response),smootho2));
                    delta_integ(ori,sf,tf) = trapz(smooth(mean_plotdata{ori,sf,tf}(response),smootho2));
                    delta_meanF0(ori,sf,tf) = nanmean(mean_plotdata{ori,sf,tf}(baseline));
                    delta_sigmaF0(ori,sf,tf) = nanstd(mean_plotdata{ori,sf,tf}(baseline));
                    delta_meanrate(ori,sf,tf) = sum(sum_plotdata{ori,sf,tf}(response)) / reps / stimdur;
                    
                    if plotdata
                        
                        hline(0, ':k')
                        
                        %send patch to back
                        h = get(gca,'Children');
                        p = find(h==findobj(h, 'Tag','Stim'));
                        
                        if ~isempty(p)
                            set(gca,'Children',[h(h~=findobj(h, 'Tag','Stim')); h(p)]);
                            
                            if rec == 1
                                plotpos = get(gca, 'Position');
                                %                         text(diff(x)/2, y(1)-0.1*y(1), x_sym_strings(ori), 'FontName', 'OriSymbols', 'FontSize',20,'LineStyle', 'none', 'Color', x_col(ori));
                                if size(x_col,1)==1
                                    text(diff(x)/2, y(1)-0*y(1), x_sym_strings(ori), 'FontName', 'OriSymbols', 'FontSize',35,'LineStyle', 'none', 'Color', x_col(ori));
                                else
                                    text(diff(x)/2, y(1)-0*y(1), x_sym_strings(ori), 'FontName', 'OriSymbols', 'FontSize',35,'LineStyle', 'none', 'Color', x_col(ori,:));
                                end
                            end
                            
                            
                        end
                        %                  text(diff(x)/2, y(1)-0.1*y(1), roi_data{rec}.fnames, 'FontName', 'OriSymbols', 'FontSize',20,'LineStyle', 'none', 'Color', x_col(ori,:));
                        if ori == oris && ey == eyes
                            
                            %draw cells
                            if ~isempty(ROIimgs_ovl)
                                dispovl = ROIimgs_ovl{rec}.ROI{ROI(rec)};
                                dispodi = ROIimgs_odi{rec}.ROI{ROI(rec)};
                                subplot_q(recs,oris*eyes + furtherfig,k);
                                imagesc(dispovl); truesize ; colormap(gray)
                                set(gca, 'visible', 'off');
                                set(gca, 'Color', 'none');
                                subplot_q(recs,oris*eyes + furtherfig,k+1);
                                imagesc(dispodi); truesize
                                set(gca, 'visible', 'off');
                                set(gca, 'Color', 'none');
                            end
                            
                            %                             k = k + furtherfig;
                        end
                        
                    end
                end
            end
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
            
            yl = ylim;
            y = [yl(2) yl(2)];
            y2 = [yl(1) yl(1)];
            p = findobj(gcf, 'Tag','Stim');
            set(p(:), 'ydata', [y fliplr(y2)]);
            
            p2 = findobj(gcf, 'Tag','scalebarx');
            set(p2(:), 'ydata', [[miny miny] [miny+tyl(2)* 0.01 miny+tyl(2)* 0.01]]);
            
            t = findobj(gcf, 'Type','text');
            poss = get(t(:), 'Position');
            for u = 1:length(poss);
                set(t(u), 'position', [poss{u}(1)+0.5  tyl(2) - 0.1 * tyl(2) 0],'FontSize',35)
            end
        end
    end
end

if plotdata
    tightfig;
    movegui('center');
end

roinum = roinum + 1;
% fillPage(gcf, 'papersize', 'A4'); orient landscape;
