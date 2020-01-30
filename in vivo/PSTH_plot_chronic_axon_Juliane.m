function [DFpeak F02  SigF0 delta_max delta_mean delta_integ delta_meanF0 delta_sigmaF0 mean_plotdata plot_data_all delta_meanrate smoothrate b a nerds int] = PSTH_plot_chronic(roi_data, ana_data, ids, ROI, ratio, single_roi, plotdata, npfct, eyes, basel, runexclude, varargin)

% disp(['extracting PSTH data of rec: ' ana_data{1}.info.FileID])
% plotdata = 0;
global roinum
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
% eyes = 2;
if isempty(npfct)
    npfct = 0.7;
end

if isempty(varargin)
    ROIimgs_ovl = []; %!!!ONLY RUN FOR VERY FIRST ANALYSIS OF ROIS!
    %     ROIimgs_green = [];
    ROIimgs_odi = [];
    
    axondf_already = 0;
else
    ROIimgs_ovl = varargin{1};
    %     ROIimgs_green = varargin{2};
    ROIimgs_odi = varargin{2};
    
    axondf_already = 0;
end

npfct_extract = 0;

pubfig = 1;
level = 4;

% plotdata = 1;

np_subtract = 0;

if ~np_subtract || npfct==0
    disp('No Neuropil Subtraction!')
end

poststim = 5;
prestim = 3;
stimdur = 5;
stimcutoff = 0.1;
F0per = 1;

episodic_F0 = 1;
detrend = 0;
useF0mode = 0;
useF0mininterp = 1;

df2nd           = 1;     % trace-by-trace secondary DF correction
noneg           = 1;     % cut off at 0 ;

ylimover = 0;
ylo = [-200 1000];

% spike inference
smc_oopsi = 0;
nerdsrun = 0;

% hill correction
hillcorrect = 0;

% NP parameters
medianadd = 1;

endper = 1;
endsub = 0;
startper = .5;
timewidth = 0.5; %timewidth of filter in seconds

% bleedthrough corr
bleedcorr = 0;
disp('NOT correcting for g->R bleedthrough');
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

if smc_oopsi
    miny = 0;
    %     plotdata = 0;
else
    a = [];
    b = [];
end
minx = -1;

scalebarx = 10;
scalebary = 200;

if smc_oopsi
    scalebary = 0.3;
end

stims = size(unique(ids{1}(1).stimseq_deg'),1);

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
    %
    %     pos = get(0, 'ScreenSize');
    %     pos(4) = (pos(4)-111)/2;
    %     pos(2) = 45+pos(4);
    %     %     set(gcf, 'Position' ,pos,  'Color', 'n');
    %     set(gcf, 'Position' ,pos);
    
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
    
    %
end


for rec = 1:recs;
    
    SamplingFreq1 = regexp(ana_data{rec}.info.ImageDescription, '(?<=scanFrameRate = )\d+\.?\d*', 'match');
    SamplingFreq1 = str2num(SamplingFreq1{1}) / level;
    
    smootho = round(timewidth * SamplingFreq1);
    cutoff = (0.44294/(sqrt(smootho^2-1)))*SamplingFreq1;
    %     disp(['moving average 3db cutoff: ' num2str(cutoff) ' Hz']);
    
    response = round(prestim*SamplingFreq1):round((prestim+stimdur-stimcutoff)*SamplingFreq1);
    baseline = round((prestim-F0per)*SamplingFreq1):round((prestim)*SamplingFreq1);
    endline = round((endper)*SamplingFreq1);
    startline = round((startper)*SamplingFreq1);
    
    if ~axondf_already
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
        %         if noneg;
        %             roi_data{rec}.ROIs(ROI(rec)).activity(roi_data{rec}.ROIs(ROI(rec)).activity <=0) = 0.01;
        %             roi_data{rec}.ROIs(ROI(rec)).activity_r(roi_data{rec}.ROIs(ROI(rec)).activity_r <=0) = 0.01;
        %             roi_data{rec}.np(ROI(rec)).activity(roi_data{rec}.np(ROI(rec)).activity <=0) = 0.01;
        %             roi_data{rec}.np(ROI(rec)).activity_r (roi_data{rec}.np(ROI(rec)).activity_r <=0) = 0.01;
        %         end
        if bleedcorr
            
            roi_data{rec}.ROIs(rec).activity_r = roi_data{rec}.ROIs(rec).activity_r - green2red_bleed * roi_data{rec}.ROIs(rec).activity;
            roi_data{rec}.np(rec).activity_r = roi_data{rec}.np(rec).activity_r - green2red_bleed * roi_data{rec}.np(rec).activity;
        end
        
    end
    
    
    
    for ey = 1:eyes;
        
        ROIs = length(roi_data{rec}.ROIs);
        oris = size(ids{rec}(ey).stim_boundaries,1);
        reps = size(ids{rec}(ey).stim_boundaries,2);
        
        plot_data = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));
        
        raw_data = roi_data{rec}.ROIs(ROI(rec)).activity(ids{rec}(ey).eye_open);
        full_raw_data = roi_data{rec}.ROIs(ROI(rec)).activity;
        
        try
            np_data = roi_data{rec}.np(ROI(rec)).activity(ids{rec}(ey).eye_open);
        catch
            %             disp('NO NEUROPIL! Using flatline!');
            roi_data{rec}.np(ROI(rec)).activity = zeros(size(full_raw_data));
            np_data = roi_data{rec}.np(ROI(rec)).activity(ids{rec}(ey).eye_open);
        end
        
        if ~np_subtract
            roi_data{rec}.np(ROI(rec)).activity = zeros(size(full_raw_data));
            np_data = roi_data{rec}.np(ROI(rec)).activity(ids{rec}(ey).eye_open);
        end
        
        full_np_data = roi_data{rec}.np(ROI(rec)).activity;
        
        if ratio
            plot_data_r = NaN(reps,ceil((prestim+stimdur+poststim)*SamplingFreq1));
        end
        
        if ratio
            raw_data_r = roi_data{rec}.ROIs(ROI(rec)).activity_r(ids{rec}(ey).eye_open);
            np_data_r = roi_data{rec}.np(ROI(rec)).activity_r(ids{rec}(ey).eye_open);
            
            full_raw_data_r = roi_data{rec}.ROIs(ROI(rec)).activity_r;
            full_np_data_r = roi_data{rec}.np(ROI(rec)).activity_r;
            %                     [~, prctrace_r] = psmooth(raw_data_r(100:end-500));
            %                     [~, prctrace_rnp] = psmooth(np_data_r(100:end-500));
        end
        
        % neuropil correction (including median of np_trace
        % addition!)
        
        if npfct_extract
            %estimating the linear neuropil contamination factor similar to
            %the dendrite contribution in spine imaging data using
            %robustfit on the NP vs. cell samplewise plot.
            
            [int] = robustfit([np_data(200:end)], [raw_data(200:end)], 'andrews', [], 'off' );
            %              [int] = robustfit([np_data(200:end)], [raw_data(200:end)], 'andrews', [], 'on' );
            
            %             npfct = int(2);
            npfct = int;
            %             inter = int(1);
            inter = 0;
            
            %             figure(21414)
            %             plot(np_data(100:end), raw_data(100:end), '.k'); hold on;
            %             plot(np_data(100:end), inter + npfct * np_data(100:end),'r','LineWidth',2);
            %             xlabel('NP signal raw (F)')
            %             ylabel('Cell signal raw (F)')
            disp(['slope ' num2str(npfct) 'intercept ' num2str(inter)]);
            
        else
            inter = 0;
        end
        if ~axondf_already
            if medianadd
                raw_data = ( raw_data - npfct * np_data) + npfct * nanmedian(np_data);
                full_raw_data = ( full_raw_data - npfct * full_np_data) + npfct * nanmedian(full_np_data);
            else
                raw_data = ( raw_data - npfct * np_data);
                full_raw_data = ( full_raw_data - npfct * full_np_data);
            end
            
            if ratio
                if medianadd
                    raw_data_r = (raw_data_r - npfct * np_data_r) + npfct * nanmedian(np_data_r);
                    full_raw_data_r = (full_raw_data_r - npfct * full_np_data_r) + npfct * nanmedian(full_np_data_r);
                else
                    raw_data_r = (raw_data_r - npfct * np_data_r);
                    full_raw_data_r = (full_raw_data_r - npfct * full_np_data_r);
                end
                % ratio
                raw_data = raw_data ./ raw_data_r;
                full_raw_data = full_raw_data ./ full_raw_data_r;
            end
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
        
        if ~axondf_already
            for ori = 1:oris;
                for rep = 1:reps
                    F0_new(ori,rep) = nanmedian(raw_data(ids{rec}(ey).stim_boundaries(ori,rep,2)-1-round(F0per*SamplingFreq1):ids{rec}(ey).stim_boundaries(ori,rep,2)-1)); %TR14 -1! before the first stimframe was included - that's of course wrong!
                    F0 = nanmedian(F0_new(:));
                end
            end
        else
            F0_new = zeros(oris,reps);
            F0 = 0;
        end
        if hillcorrect
            if detrend
                [raw_data, ~] = psmooth(raw_data);
            end
            
            Fsat = 700 * 100/85;
            P.n     = 2.9 %Hill coefficient
            P.k_d   = 144  %dissociation constant
            
            theta_F = ( (raw_data - F0)/F0 * 100) / Fsat;
            theta_F_full = ( (full_raw_data - F0)/F0 * 100) / Fsat;
            
            raw_data = invHill_v2([],theta_F);
            full_raw_data = invHill_v2([],theta_F_full);
        end
        
        
        if smc_oopsi && ey ==1; %spike-rate inference based on Vogelstein's stuff
            %                     deltatrace = (raw_data - F0) ./ F0;
            
            %             full_raw_data = raw_data
            
            %remove NaN
            full_raw_data_sp = full_raw_data;
            nani = find(isnan(full_raw_data))
            
            if ~isempty(nani);
                %delete NaN
                full_raw_data_sp(nani) = [];
            end
            
            
            V.fast_iter_max = 50;
            V.smc_iter_max = 0;
            V.preprocess = 1;
            V.dt = 1/SamplingFreq1;
            V.T = size(full_raw_data_sp,2);
            V.plot = 0;
            V.fast_plot = 0;
            V.save2  = 1;
            if V.save2
                mkdir('oopsi')
            end
            V.name = ['\oopsi\ROI' num2str(roinum) 'Group' num2str(roi_data{1}.ROIs.group)];
            spikesori = 20;
            spikes = spikesori*oris*reps;
            
            P.k = log(-log(1-spikes/V.T)/V.dt);
            P.A     = 10 %1AP deltaF/F0
            P.n     = 2.9 %Hill coefficient
            P.k_d   = 144  %dissociation constant
            P.C_0   = 0;
            
            
            if  ~V.smc_iter_max == 0;
                [a b c] = run_oopsi(full_raw_data_sp, V, P);
            else
                [a b ] = run_oopsi(full_raw_data_sp, V, P);
                b.E.nbar = zeros(size(full_raw_data_sp));
            end
            
            
            if nerdsrun
                
                runs = 10;
                
                opts = NerdsPrepare(SamplingFreq1, runs)
                [gen_atom_out, spike_idx, x_hat_out, e_hat_out] = compute_nerds(full_raw_data, opts);
                recfluo = conv(gen_atom_out, x_hat_out);
                recfluo = recfluo(1:length(full_raw_data));
                raw_spikes_ev_nerds = zeros(size(recfluo));
                raw_spikes_ev_nerds(spike_idx{1}) = 1;
                raw_spikes_ev_nerds  =logical(raw_spikes_ev_nerds);
            end
            
            if ~isempty(nani);
                %re-add previous NaN with 0
                full_raw_data_sp = [full_raw_data_sp(1:nani(1)-1) zeros(1,length(nani)) full_raw_data_sp(nani(1):end)];
                b.E.nbar = [b.E.nbar(1:nani(1)-1) NaN(1,length(nani)) b.E.nbar(nani(1):end)];
                a.n = [a.n(1:nani(1)-1); NaN(length(nani),1); a.n(nani(1):end)];
                
                if nerdsrun
                    gen_atom_out = [gen_atom_out(1:nani(1)-1) zeros(1,length(nani)) gen_atom_out(nani(1):end)];
                    x_hat_out= [x_hat_out(1:nani(1)-1) zeros(1,length(nani)) x_hat_out(nani(1):end)];
                    e_hat_out= [e_hat_out(1:nani(1)-1) zeros(1,length(nani)) e_hat_out(nani(1):end)];
                end
            end
            
            if nerdsrun
                recfluo = conv(gen_atom_out, x_hat_out);
                recfluo = recfluo(1:length(full_raw_data));
                raw_spikes_ev_nerds = zeros(size(recfluo));
                raw_spikes_ev_nerds(spike_idx{1}) = 1;
                raw_spikes_ev_nerds  =logical(raw_spikes_ev_nerds);
            end
            
            zlev = 2;
            zlev_s = 2;
            
            raw_spikes=find(b.E.nbar>zlev_s*nanstd(b.E.nbar));
            raw_spikes_foopsi=find(a.n>zlev*nanstd(a.n));
            raw_spikes_ev=b.E.nbar>zlev_s*nanstd(b.E.nbar);
            raw_spikes_ev_foopsi=a.n>zlev*nanstd(a.n);
            
            timebase = (1:length(full_raw_data_sp)) ./ SamplingFreq1;
            plotrate = 0;
            
            %convert to spikerate
            sigma = SamplingFreq1/2;
            edges  = [-3*sigma:1:3*sigma];
            kernel = normpdf(edges,0,sigma);
            
            smoothrate = conv(single(raw_spikes_ev), kernel);
            smoothrate_foopsi = conv(single(raw_spikes_ev_foopsi), kernel);
            if nerdsrun
                smoothrate_nerds = conv(single(raw_spikes_ev_nerds), kernel);
            end
            center = ceil(length(edges)/2);
            
            smoothrate = smoothrate(center:length(timebase)+center-1);
            smoothrate_foopsi = smoothrate_foopsi(center:length(timebase)+center-1);
            if nerdsrun
                smoothrate_nerds = smoothrate_nerds(center:length(timebase)+center-1);
            end
            runavrate = smooth(raw_spikes_ev,SamplingFreq1);
            runavrate_foopsi = smooth(raw_spikes_ev_foopsi,SamplingFreq1);
            
            if plotrate
                if nerdsrun
                    subf = 1;
                else
                    subf = 0;
                end
                try
                    close(2326463);kkk = figure(2326463);
                catch
                    kkk = figure(2326463);
                end
                %             h = stem(timebase(raw_spikes), ones(size(raw_spikes))*max(full_raw_data),'-', 'Color', [0.75 0.75 0.75]); set(h, 'Marker', 'none'); hold on
                subplot(1+subf,1,1)
                h = stem(timebase(raw_spikes_ev_foopsi), ones(size(find(raw_spikes_ev_foopsi==1)))*max(full_raw_data_sp),'-', 'Color', [0.5 0.5 0.5]); set(h, 'Marker', 'none'); hold on
                %                 plot(timebase,runavrate, '-b');
                %                 plot(timebase,smoothrate, '-r');
                plot(timebase,smoothrate_foopsi, '-g');
                
                plot(timebase,full_raw_data_sp, 'k'); hold on;
                title('fast_oopsi');
                %                 tightfig;
                set(kkk, 'Position', [          647         125        1311        1270]);
                
                if nerdsrun
                    subplot(1+subf,1,1+subf)
                    h = stem(timebase(raw_spikes_ev_nerds), ones(size(find(raw_spikes_ev_nerds==1)))*max(full_raw_data_sp),'-', 'Color', [0.5 0.5 0.5]); set(h, 'Marker', 'none'); hold on
                    %                 plot(timebase,runavrate, '-b');
                    %                 plot(timebase,smoothrate, '-r');
                    plot(timebase,smoothrate_nerds, '-g');
                    plot(timebase, recfluo + F0, '-r')
                    plot(timebase,full_raw_data_sp, 'k'); hold on;
                    title('NERDS');
                    
                    %consolidate output
                    nerds.gen_atom_out =gen_atom_out;
                    nerds.spike_idx = spike_idx;
                    nerds.x_hat_out = x_hat_out;
                    nerds.e_hat_out = e_hat_out;
                    nerds.F0 = F0;
                    nerds.recfluo = recfluo
                    xlabel(['ROI' num2str(roinum) 'Group' num2str(roi_data{1}.ROIs.group) ' time [s]']);
                    set(gcf, 'name', ['ROI' num2str(roinum) 'Group' num2str(roi_data{1}.ROIs.group) ]);
                    saveas(gcf, [cd '\' V.name '.tif']);
                    saveas(gcf, [cd '\' V.name '.fig']);
                    %                 tightfig;
                    set(kkk, 'Position', [          647         125        1311        1270]);
                end
                
                %                 figure(23553);
            end
            
        end
        
        if smc_oopsi
            if  V.smc_iter_max == 0;
                %                 raw_data = single(raw_spikes_ev_foopsi(ids{rec}(ey).eye_open))';
                raw_data = single(a.n(ids{rec}(ey).eye_open))'; % using foopsi_n
            else
                raw_data = single(raw_spikes_ev(ids{rec}(ey).eye_open))';
            end
        end
        
        if ~axondf_already
            if ~episodic_F0
                % %             detrend data
                if ~smc_oopsi && detrend
                    [raw_data, t] = psmooth(raw_data);
                end
                
                for ori = 1:oris;
                    % %             generate episodic prestim F0 data and take median of these for overall F0 value of detrended raw data
                    for rep = 1:reps
                        F0_new(ori,rep) = nanmedian(raw_data(ids{rec}(ey).stim_boundaries(ori,rep,2)-1-round(F0per*SamplingFreq1):ids{rec}(ey).stim_boundaries(ori,rep,2)-1)); %TR14 -1! before the first stimframe was included - that's of course wrong!
                        F0 = nanmedian(F0_new(:));
                        if useF0mode
                            disp('Using mode of F0 distributions as global F0');
                            F0 = mode(F0_new(:));
                        end
                        
                        
                        %                     if useconcatbase
                        %                         F0_concat_base = getF0frombaseline(ids_new)
                        %                     end
                    end
                end
                
            end
            
            if useF0mininterp
                [F0_min_new, ~,F0_min_interp,~] = getF0frombaseline(raw_data, ids{rec}, F0per, SamplingFreq1);
            end
            if noneg;
                raw_data(raw_data <=0) = 0.01;
                %             roi_data{rec}.ROIs(ROI(rec)).activity_r(roi_data{rec}.ROIs(ROI(rec)).activity_r <=0) = 0.01;
                %             roi_data{rec}.np(ROI(rec)).activity(roi_data{rec}.np(ROI(rec)).activity <=0) = 0.01;
                %             roi_data{rec}.np(ROI(rec)).activity_r (roi_data{rec}.np(ROI(rec)).activity_r <=0) = 0.01;
            end
            
            
        else
            F0_min_new = F0_new;
        end
        
        
        for ori = 1:oris;
            
            for rep = 1:reps
                
                ids{rec}(ey).stim_boundaries(:,:,4) = ids{rec}(ey).stim_boundaries(:,:, 3) + floor(SamplingFreq1*poststim);
                ids{rec}(ey).stim_boundaries(ids{rec}(ey).stim_boundaries>length(roi_data{rec}.ROIs(ROI(rec)).activity(ids{rec}(ey).eye_open))) = length(roi_data{rec}.ROIs(ROI(rec)).activity(ids{rec}(ey).eye_open));
                
                if ~axondf_already;
                    if episodic_F0 && ~useF0mininterp
                        F0 = nanmedian(raw_data(ids{rec}(ey).stim_boundaries(ori,rep,2)-1-round(F0per*SamplingFreq1):ids{rec}(ey).stim_boundaries(ori,rep,2)-1)); %TR14 -1! before the first stimframe was included - that's of course wrong!
                    elseif episodic_F0 && useF0mininterp
                        F0 =  F0_min_new(ori,rep,ey);
                    end
                else
                    F0 = 0;
                end
                
                %                     F0(F0<1)=1;
                DF = raw_data(ids{rec}(ey).stim_boundaries(ori,rep,2)-round(prestim*SamplingFreq1):ids{rec}(ey).stim_boundaries(ori,rep,4))-F0;
                
                if ~smc_oopsi
                    plot_data(rep,1:length(DF)) =(DF./F0) * 100;
                else
                    plot_data(rep,1:length(DF)) =(DF);
                    %                      plot_data(rep,1:length(DF)) =(DF./F0) * 100;
                end
                
                if runexclude && ~plotdata
                    if  ids{rec}(ey).running_maxspeed(ori,rep) > speed_crit
                        %                         ids{rec}(ey).running_dist(ori,rep) > dist_crit
                        plot_data(rep,1:length(DF)) = NaN;
                        %                          disp(['Recording ' num2str(rec) ' eye ' num2str(ey) ' ori ' num2str(ori) ' rep ' num2str(rep) ' is a running trial -> NaN!'])
                    end
                end
                
                timebase = (1:length(plot_data)) / SamplingFreq1;
                
                %% consolidate data output fdrom single trace
                DFpeak(ori,rep,ey) = nanmax(plot_data(rep,response));
                F02(ori,rep,ey) = nanmean(plot_data(rep,baseline));
                SigF0(ori,rep,ey) = nanstd(plot_data(rep,baseline));
                Fend(ori,rep,ey) = nanmean(plot_data(rep,end-endline:end));
                Fstart(ori,rep,ey) = nanmean(plot_data(rep,1:startline));
                
                %correct for corrupted baseline by adding the post-decay
                if endsub
                    if Fstart(ori,rep,ey) < - 6 * SigF0(ori,rep,ey);
                        DFpeak(ori,rep,ey) = DFpeak(ori,rep,ey) + abs(Fstart(ori,rep,ey));
                        F02(ori,rep,ey) = F02(ori,rep,ey);
                        plot_data(rep,1:length(DF)) = plot_data(rep,1:length(DF)) + abs(Fstart(ori,rep,ey));
                        disp('F0 corrupted - OFFSET: Fend!!!!!!!!');
                    end
                end
                
                if plotdata
                    hdl(ey,ori) = subplot_q(recs,oris*eyes + furtherfig,k);
                    if ~smc_oopsi
                        ptr = plot(timebase(1:downsamp:end),plot_data(rep,1:downsamp:end), 'k');
                        try
                            if  ids{rec}(ey).running_maxspeed(ori,rep) > speed_crit
                                %                             ids{rec}(ey).running_dist(ori,rep) > dist_crit
                                set(ptr, 'Color', 'b');
                            end
                        end
                        %                          disp(['Recording ' num2str(rec) ' eye ' num2str(ey) ' ori ' num2str(ori) ' rep ' num2str(rep) ' is a running trial -> NaN!'])
                        
                    else
                        ptr = plot(timebase(1:downsamp:end),plot_data(rep,1:downsamp:end), 'Color', [0 0 0]);
                        
                    end
                    hold on;
                end
            end
            
            if runexclude
                usetrials = find(~isnan(DFpeak(ori,:,ey)));
                if length(usetrials) < mintrials
                    %                     disp(['Recording ' num2str(rec) ' eye ' num2str(ey) ' ori ' num2str(ori) 'has more than ' num2str(mintrials) ' running trials -> Exclude timepoint!'])
                    %                   usetrials = [usetrials randsample(usetrials,mintrials - length(usetrials))];
                    % take all trials then and live with the NaNs (sort them
                    % out later!_
                    usetrials = randsample(length(DFpeak(ori,:,ey)),mintrials);
                else
                    usetrials =  randsample(usetrials,mintrials);
                end
                excludetrials =  setdiff(1:reps, usetrials);
                
                DFpeak_re(ori,:,ey) =  DFpeak(ori,usetrials,ey);
                F02_re(ori,:,ey) =  F02(ori,usetrials,ey);
                SigF0_re(ori,:,ey) =  SigF0(ori,usetrials,ey);
                Fend_re(ori,:,ey) =  Fend(ori,usetrials,ey);
                Fstart_re(ori,:,ey) =  Fstart(ori,usetrials,ey);
            else
                usetrials = 1:length(DFpeak(ori,:,ey));
                
                
            end
            
            
            %             mean_plotdata{ori,ey} = smooth(nanmean(plot_data(usetrials,:),1), smootho);
            mean_plotdata{ori,ey} = nanmean(plot_data(usetrials,:),1);
            plot_data_all{ori,ey} = plot_data;
            sum_plotdata{ori,ey} = sum(plot_data,1);
            
            timebase = (1:length(mean_plotdata{ori,ey})) / SamplingFreq1;
            
            if plotdata
                if ~smc_oopsi
                    plot(timebase(1:downsamp:end), mean_plotdata{ori,ey}(1:downsamp:end), 'r', 'LineWidth', 1);
                else
                    %                     [histdat thist] = hist(find(sum_plotdata{ori,ey}(1:downsamp:end)),[1:SamplingFreq1:length(sum_plotdata{ori,ey}(1:downsamp:end))]);
                    %                     stairs(thist ./ SamplingFreq1,histdat,  'r', 'LineWidth', 1);
                    plot(timebase(1:downsamp:end), mean_plotdata{ori,ey}(1:downsamp:end), 'r', 'LineWidth', 1);
                end
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
                    if smc_oopsi
                        lazy = [0 0.3];
                    end
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
                %             axis tight;
                
                yl = ylim;
                x = [prestim prestim+stimdur];
                y = [yl(2) yl(2)];
                y2 = [yl(1) yl(1)];
                
                if ey ==1 ;
                    patch([x fliplr(x)],[y fliplr(y2)], 'r','facealpha',0.1,...
                        'edgecolor','none', 'Tag', 'Stim');
                elseif ey == 2 ;
                    patch([x fliplr(x)],[y fliplr(y2)], 'b','facealpha',0.1,...
                        'edgecolor','none', 'Tag', 'Stim');
                elseif ey == 3;
                    patch([x fliplr(x)],[y fliplr(y2)], 'w','facealpha',0.1,...
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
                    
                    %MD marker
                    if ~isempty(basel)
                        if rec == basel + 1
                            set(ptsb1, 'FaceColor', 'r');
                            set(ptsb2, 'FaceColor', 'r');
                        end
                    end
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
            
            delta_max(ori,ey) = nanmax(smooth(mean_plotdata{ori,ey}(response),smootho2));
            delta_mean(ori,ey) = nanmean(smooth(mean_plotdata{ori,ey}(response),smootho2));
            delta_integ(ori,ey) = trapz(smooth(mean_plotdata{ori,ey}(response),smootho2));
            delta_meanF0(ori,ey) = nanmean(mean_plotdata{ori,ey}(baseline));
            delta_sigmaF0(ori,ey) = nanstd(mean_plotdata{ori,ey}(baseline));
            delta_meanrate(ori,ey) = sum(sum_plotdata{ori,ey}(response)) / reps / stimdur;
            
            if plotdata
                if ~pubfig && ~ smc_oopsi
                    hline(delta_max(ori, ey), '-r');
                end
                if smc_oopsi;
                    %                     hline(delta_meanrate(ori, ey), '-r');
                    hline(delta_max(ori, ey), '-r');
                end
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
                            text(diff(x)/2, y(1)-0*y(1), x_sym_strings(ori), 'FontName', 'OriSymbols', 'FontSize',30,'LineStyle', 'none', 'Color', x_col(ori));
                        else
                            text(diff(x)/2, y(1)-0*y(1), x_sym_strings(ori), 'FontName', 'OriSymbols', 'FontSize',30,'LineStyle', 'none', 'Color', x_col(ori,:));
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
                    if  isfield(ana_data{rec}, 'Fit')
                        tunq(rec) = subplot_q(recs,oris*eyes + furtherfig,k+2:k+3);
                        hfitc(rec) = plot([ana_data{rec}.Fit(ROI(rec)).contra.FittedData],'b', 'LineWidth', 1); hold on;
                        hfiti(rec) = plot([ana_data{rec}.Fit(ROI(rec)).ipsi.FittedData],'r', 'LineWidth', 1);
                        if eyes == 3;
                            hfitb(rec) = plot([ana_data{rec}.Fit(ROI(rec)).bino.FittedData],'k', 'LineWidth', 1);
                        end
                        
                        oppdir_contra = mod(round([ana_data{rec}.Fit(ROI(rec)).contra.PrefDir]+(length([ana_data{rec}.Fit(ROI(rec)).contra.FittedData])/2))-1, length([ana_data{rec}.Fit(ROI(rec)).contra.FittedData]))+1;
                        oppdir_ipsi = mod(round([ana_data{rec}.Fit(ROI(rec)).ipsi.PrefDir]+(length([ana_data{rec}.Fit(ROI(rec)).ipsi.FittedData])/2))-1, length([ana_data{rec}.Fit(ROI(rec)).ipsi.FittedData]))+1;
                        if eyes == 3;
                            oppdir_bino = mod(round([ana_data{rec}.Fit(ROI(rec)).bino.PrefDir]+(length([ana_data{rec}.Fit(ROI(rec)).bino.FittedData])/2))-1, length([ana_data{rec}.Fit(ROI(rec)).bino.FittedData]))+1;
                        end
                        
                        hfitcp(rec) = plot([ana_data{rec}.peaks(ROI(rec)).oris], [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_contra], 'ob', 'MarkerFaceColor', 'w');
                        hfitip(rec) = plot([ana_data{rec}.peaks(ROI(rec)).oris], [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_ipsi], 'or','MarkerFaceColor', 'w');
                        if eyes == 3;
                            hfitbp(rec) = plot([ana_data{rec}.peaks(ROI(rec)).oris], [ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_bino], 'ok','MarkerFaceColor', 'w');
                        end
                        
                        scaleamp(rec,1) = max([ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_contra]);
                        scaleamp(rec,2) = max([ana_data{rec}.peaks(ROI(rec)).deltapeaks_averagetrace_ipsi]);
                        
                        %                             lazy = [-20 max(scaleamp(:)) + 0 * max(scaleamp(:))];
                        %                             lazy = [miny nanmax(displaypeaks(:))+0.1*nanmax(displaypeaks(:))];
                        set(tunq(rec), 'ylim', lazy, 'xlim', [0 length([ana_data{rec}.Fit(ROI(rec)).contra.FittedData])], 'XTick', [ana_data{rec}.peaks(ROI(rec)).oris 360])
                        set(gca, 'visible', 'on');
                        
                        l3(rec)= vline([ana_data{rec}.Fit(ROI(rec)).contra.PrefDir], ':b', num2str(round([ana_data{rec}.Fit(ROI(rec)).contra.PrefDir])));
                        l4(rec) = vline([ana_data{rec}.Fit(ROI(rec)).ipsi.PrefDir], ':r', num2str(round([ana_data{rec}.Fit(ROI(rec)).ipsi.PrefDir])));
                        l6(rec) = vline(oppdir_contra, ':b',num2str(round(oppdir_contra)));
                        l7(rec)=vline(oppdir_ipsi, ':r', num2str(round(oppdir_ipsi)));
                        
                        if eyes == 3;
                            l5(rec) = vline([ana_data{rec}.Fit(ROI(rec)).bino.PrefDir], ':k', num2str([ana_data{rec}.Fit(ROI(rec)).bino.PrefDir]));
                        end
                        
                        hline(0, ':k')
                    end
                    k = k + furtherfig;
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
            if~smc_oopsi
                set(hdl(:), 'ylim', [miny nanmax(delta_max(:)) + 0.333*nanmax(delta_max(:))]);
                if ylimover
                    set(hdl(:), 'ylim', ylo);
                end
            else
                %                 set(hdl(:), 'ylim', [miny nanmax(delta_meanrate(:)) + 0.333*nanmax(delta_meanrate(:))]);
                set(hdl(:), 'ylim', [miny nanmax(delta_max(:)) + 0.333*nanmax(delta_max(:))]);
                if ylimover
                    set(hdl(:), 'ylim', ylo);
                end
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
                set(t(u), 'position', [poss{u}(1)+0.5  tyl(2) - 0.1 * tyl(2) 0],'FontSize',40)
            end
        end
    end
end

if plotdata
    tightfig;
end

if ~smc_oopsi;
    delta_meanrate =[];
    smoothrate =[];
    smc_oopsi = [];
    f_oopsi = [];
    nerds = [];
    a=[];
    b = [];
    
end

if ~nerdsrun
    nerds = [];
end
roinum = roinum + 1;
