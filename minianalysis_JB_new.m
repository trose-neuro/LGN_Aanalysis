function  [red_failure_AMPA, blue_failure_AMPA, red_failure_NMDA, blue_failure_NMDA] = ...
    minianalysis_JB_new(list, idx, clamp, pathName, fc, show, raw_traces, user, filterephys)
%Function that extracts minis by using the std threshold criterion

%list=      information of folders for each cell SW000XX
%idx=       vector with idx of which recording is mini recording
%pathName=  folder name of cell
%fc=        factor of how many stds the signal should be included
%show=      show plots or not (1 or 0)

base_start          =   1;
base_end            =   99;
pulse_start         =   100;
pulse_end           =   110;
pulse2_start        =   351;
pulse2_end          =   360;

%Calibration curves for blue/red at the 2 setups (in mW/mm2) irradiance
%compare to Klapoetke 2014
%SW
% yirr_red=(12.19*PD1-0.4319)/100;
% yirr_blue=(7.232*PD2-0.9951)/100;
% %MF
% yirr_red=(104.1 *PD1-3.467)/100;
% yirr_blue=(679.2*PD2-26.82)/100;

%% TR2019: filtering
% filterephys = 1;        % filtering yes/no? (now as function input!)
cutoff      = 500;      % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel. Use Bessel at > 4 order to prevent ripples)


%% TR2019: plot specs
plotlength = .3; %seconds
savefig = 0; %save main figure

if filterephys
    disp('- - - - - - - -')
    disp(['Filtering: ' num2str(order) ' pole ' type '-Filter w/ ' num2str(cutoff) ' Hz cutoff']);
    disp('- - - - - - - -')
end

vclamp=clamp(idx);

for i=1:length(idx)
    load([char(pathName) filesep list(idx(i)).name],'-mat');
    sr = header.ephys.ephys.sampleRate; %check sample rate
    srF = 1/(1000/sr);
    samples_per_sweep = header.ephys.ephys.traceLength*sr;
    timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
    
    ephystraces=data.ephys.trace_1;
    
%     if filterephys % TR2019: filtering
%         ephystraces = lowpassfilt(ephystraces, order, cutoff, sr, type);
%     end
    
    traces=reshape(ephystraces, samples_per_sweep, length(ephystraces)/samples_per_sweep);
    photodiode=data.acquirer.trace_1;
    photodiode=reshape(photodiode, samples_per_sweep, length(ephystraces)/samples_per_sweep);

    bs=traces(base_start*srF:base_end*srF,:);
    bs_std=std(bs);
    bs_traces=bsxfun(@minus, traces, mean(bs));
    bs_photodiode=bsxfun(@minus, photodiode, mean(photodiode(base_start*srF:base_end*srF,:)));
        
    traces_all(:,:,i)=bs_traces;
   
    double_pulse(:,i)=mean(bs_photodiode((pulse2_start+4)*srF:pulse2_end*srF))>0.025;
    
    % shape and variation of photodiode signal needs to be saved for quality control!
    temp = [mean(bs_photodiode'); std(bs_photodiode')]; 
    bs_photodiode_mean_std{i} = temp(:,1:10:end); % downsample by 10
    
    % photodiode pulse amplitude
    PD_pulse1(:,i)=mean(bs_photodiode(pulse_start *srF:pulse_end*srF,:));
    PD_pulse2(:,i)=mean(bs_photodiode(pulse2_start *srF:pulse2_end*srF,:));
    
    % Calculate irradiance from photodiode signal
    if user==0%SW 
        IR1_r(:,i)=(12.19*PD_pulse1(:,i)-0.4319)/100; 
        IR1_b(:,i)=(7.232*PD_pulse1(:,i)-0.9951)/100;
        IR2_b(:,i)=(7.232*PD_pulse2(:,i)-0.9951)/100;
    else %MF
        IR1_r(:,i)=(104.1*PD_pulse1(:,i)-3.467)/100;
        IR1_b(:,i)=(679.2*PD_pulse1(:,i)-26.82)/100;
        IR2_b(:,i)=(679.2*PD_pulse2(:,i)-26.82)/100;
    end
    
    % Test pulse for Rs,Rm and Cm calculation SW     
    if user==0
        [Rs_cell(:,i), Rm_cell(:,i), Cm_cell(:,i)]=calc_cell_acess_prop(header, bs_traces, sr, srF, base_end);
    else
        Rs_cell(:,i)=NaN;
        Rm_cell(:,i)=NaN;
        Cm_cell(:,i)=NaN;
    end
end

idx_red70=find(vclamp==1 & double_pulse==1);% Red AMPA
idx_blue70=find(vclamp==1 & double_pulse==0);% Blue AMPA
idx_red40=find(vclamp==0 & double_pulse==1);% Red NMDA
idx_blue40=find(vclamp==0 & double_pulse==0);% Blue NMDA

response_window = pulse_start*srF:pulse_end*srF+40*srF; % not sure if pulse end +40 is a good idea here
baseline_window = (base_end - (pulse_end-pulse_start+40))*srF:base_end*srF;

% Red AMPA
if ~isempty(idx_red70)
    % Add baseline subtracted traces to output vars
    if raw_traces
        red_failure_AMPA.traces_all = traces_all(1:10:end,:,idx_red70); % downsample this by 10
    else
        red_failure_AMPA.traces_all = [];
    end
    
    % add response measurements to output var
    [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = ...
        extract_events(-1*traces_all(:,:,idx_red70), baseline_window, response_window, fc, show); % trace is flipped to positive for event detection
    red_failure_AMPA.peaks = peaks; % peak amplitudes during response window
    red_failure_AMPA.test = test; % noparametric test of significance for presence of responses
    red_failure_AMPA.resp_thresh = resp_thresh; % amplitude threshold for detecting events
    red_failure_AMPA.resp_idx = resp_idx; % index of where events were detected
    red_failure_AMPA.steady_state = steady_state; % repetition at which amplitude decay has reached steadystate (estimated by fitting an exponential to event amplitudes)
    red_failure_AMPA.resp_prop = resp_prop; % response probability after steady state is reached
    red_failure_AMPA.avg_amp = avg_amp; % avg response amp after steady state is reached
   
    % add irradiance to output var
    red_failure_AMPA.photodiode_mean_std = bs_photodiode_mean_std{idx_red70}; % photodiode signal trace
    red_failure_AMPA.IR_pulse1 = IR1_r(:,idx_red70); % irradiance of first pulse
    red_failure_AMPA.IR_pulse2 = IR2_b(:,idx_red70); % irradiance of second pulse
    
    % add series resistance to output var
    red_failure_AMPA.Rs_cell = Rs_cell(:,idx_red70); % series resistance
    red_failure_AMPA.Rm_cell = Rm_cell(:,idx_red70); % membrane resistance
    red_failure_AMPA.Cm_cell = Cm_cell(:,idx_red70); % membrane capasitance
else
    red_failure_AMPA.traces_all = [];
    red_failure_AMPA.peaks = [];
    red_failure_AMPA.test = [];
    red_failure_AMPA.resp_thresh = [];
    red_failure_AMPA.resp_idx = [];
    red_failure_AMPA.steady_state = [];
    red_failure_AMPA.resp_prop = [];
    red_failure_AMPA.avg_amp = [];
    red_failure_AMPA.photodiode_mean_std = [];
    red_failure_AMPA.IR_pulse1 = [];
    red_failure_AMPA.IR_pulse2 = [];
    red_failure_AMPA.Rs_cell = [];
    red_failure_AMPA.Rm_cell = [];
    red_failure_AMPA.Cm_cell = [];
end


if ~isempty(idx_blue70)% blue AMPA
    if raw_traces
        blue_failure_AMPA.traces_all = traces_all(1:10:end,:,idx_blue70);
    else
        blue_failure_AMPA.traces_all = [];
    end
    
    [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = ...
        extract_events(-1*traces_all(:,:,idx_blue70), baseline_window, response_window, fc, show);
    blue_failure_AMPA.peaks = peaks;
    blue_failure_AMPA.test = test;
    blue_failure_AMPA.resp_thresh = resp_thresh;
    blue_failure_AMPA.resp_idx = resp_idx;
    blue_failure_AMPA.steady_state = steady_state;
    blue_failure_AMPA.resp_prop = resp_prop;
    blue_failure_AMPA.avg_amp = avg_amp; 

    blue_failure_AMPA.photodiode_mean_std = bs_photodiode_mean_std{idx_blue70};
    blue_failure_AMPA.IR_pulse1 = IR1_r(:,idx_blue70);
    
    blue_failure_AMPA.Rs_cell = Rs_cell(:,idx_blue70);
    blue_failure_AMPA.Rm_cell = Rm_cell(:,idx_blue70);
    blue_failure_AMPA.Cm_cell = Cm_cell(:,idx_blue70);
else
    blue_failure_AMPA.traces_all = [];
    blue_failure_AMPA.peaks = [];
    blue_failure_AMPA.test = [];
    blue_failure_AMPA.resp_thresh = [];
    blue_failure_AMPA.resp_idx = [];
    blue_failure_AMPA.steady_state = [];
    blue_failure_AMPA.resp_prop = [];
    blue_failure_AMPA.avg_amp = []; 
    blue_failure_AMPA.photodiode_mean_std = [];
    blue_failure_AMPA.IR_pulse1 = [];
    blue_failure_AMPA.Rm_cell = [];
    blue_failure_AMPA.Cm_cell = [];
end
    
if ~isempty(idx_red40)% Red NMDA
    if raw_traces
        red_failure_NMDA.traces_all = traces_all(1:10:end,:,idx_red40);
    else
        red_failure_NMDA.traces_all = [];
    end
    
    [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = ...
        extract_events(traces_all(:,:,idx_red40), baseline_window, response_window, fc, show);
    red_failure_NMDA.peaks = peaks;
    red_failure_NMDA.test = test;
    red_failure_NMDA.resp_thresh = resp_thresh;
    red_failure_NMDA.resp_idx = resp_idx;
    red_failure_NMDA.steady_state = steady_state;
    red_failure_NMDA.resp_prop = resp_prop;
    red_failure_NMDA.avg_amp = avg_amp; 

    red_failure_NMDA.photodiode_mean_std = bs_photodiode_mean_std{idx_red40};
    red_failure_NMDA.IR_pulse1 = IR1_r(:,idx_red40);
    red_failure_NMDA.IR_pulse2 = IR2_b(:,idx_red40);

    red_failure_NMDA.Rs_cell = Rs_cell(:,idx_red40);
    red_failure_NMDA.Rm_cell = Rm_cell(:,idx_red40);
    red_failure_NMDA.Cm_cell = Cm_cell(:,idx_red40);
    
else
    red_failure_NMDA.traces_all = [];
    red_failure_NMDA.peaks = [];
    red_failure_NMDA.test = [];
    red_failure_NMDA.resp_thresh = [];
    red_failure_NMDA.resp_idx = [];
    red_failure_NMDA.steady_state = [];
    red_failure_NMDA.resp_prop = [];
    red_failure_NMDA.avg_amp = []; 
    red_failure_NMDA.photodiode_mean_std = [];
    red_failure_NMDA.IR_pulse1 = [];
    red_failure_NMDA.IR_pulse2 = [];
    red_failure_NMDA.Rs_cell = [];
    red_failure_NMDA.Rm_cell = [];
    red_failure_NMDA.Cm_cell = [];
end

if ~isempty(idx_blue40) %blue NMDA
    if raw_traces
        blue_failure_NMDA.traces_all = traces_all(1:10:end,:,idx_blue40);
    else
        blue_failure_NMDA.traces_all = [];
    end
    
    [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = ...
        extract_events(traces_all(:,:,idx_blue40), baseline_window, response_window, fc, show);
    blue_failure_NMDA.peaks = peaks;
    blue_failure_NMDA.test = test;
    blue_failure_NMDA.resp_thresh = resp_thresh;
    blue_failure_NMDA.resp_idx = resp_idx;
    blue_failure_NMDA.steady_state = steady_state;
    blue_failure_NMDA.resp_prop = resp_prop;
    blue_failure_NMDA.avg_amp = avg_amp; 
    
    blue_failure_NMDA.photodiode_mean_std = bs_photodiode_mean_std{idx_blue40};
    blue_failure_NMDA.IR_pulse1 = IR1_r(:,idx_blue40);
    
    blue_failure_NMDA.Rs_cell = Rs_cell(:,idx_blue40);
    blue_failure_NMDA.Rm_cell = Rm_cell(:,idx_blue40);
    blue_failure_NMDA.Cm_cell = Cm_cell(:,idx_blue40);
else
    blue_failure_NMDA.traces_all = [];
    blue_failure_NMDA.peaks = [];
    blue_failure_NMDA.test = [];
    blue_failure_NMDA.resp_thresh = [];
    blue_failure_NMDA.resp_idx = [];
    blue_failure_NMDA.steady_state = [];
    blue_failure_NMDA.resp_prop = [];
    blue_failure_NMDA.avg_amp = []; 
    blue_failure_NMDA.photodiode_mean_std = [];
    blue_failure_NMDA.IR_pulse1 = [];
    blue_failure_NMDA.Rs_cell = [];
    blue_failure_NMDA.Rm_cell = [];
    blue_failure_NMDA.Cm_cell = [];
end

end

function [peaks, test, resp_thresh, resp_idx, steady_state, resp_prop, avg_amp] = extract_events(bs_traces, baseline_window, response_window, fc, show_plot)
if nargin<4
    show_plot = 0;
end

mean_resp = nanmean(bs_traces(response_window,:));
mean_base = nanmean(bs_traces(baseline_window,:));
bs_std = std(bs_traces(baseline_window,:),1);
resp_thresh = bs_std*4;

try 
    test = signrank(mean_base,mean_resp)<0.05;
catch
    test = 0;
end

peaks = max(bs_traces(response_window,:));

if test
    if show_plot
        % To see how the distribution of peak hights differse between
        % baseline and stim period
        peak_base = min(bs_traces(baseline_window,:));
        figure; hist([peaks;peak_base]',100)
    end
    
    % error occurs here, check!
    resp_idx = peaks>mean(resp_thresh); % which trials show a response 
    
    % Fit exp decay in order to only measure steady state responses
    try
        % single exponential
        sucs_peaks = peaks(find(resp_idx))';
        x = [0:sum(resp_idx)-1]';
        fo = fitoptions('Method','NonlinearLeastSquares', ...
            'Lower',[0 -1 0], ...
            'Upper', [max(sucs_peaks)*10 0 max(sucs_peaks)],...
            'StartPoint',[max(sucs_peaks)*3 -0.0005 min(sucs_peaks)]); % adjust start settings
        ft = fittype('a*exp(b*(x))+c', 'options', fo);
        [f gof]=fit(x,sucs_peaks,ft);
        yf=f.a.*exp(f.b.*x)+f.c;
        if show_plot
        	figure; plot(x,yf); hold on; plot(x,sucs_peaks,'ro')
        end
        
        R2_exp = sum((sucs_peaks-yf).^2);
        R2_linear = sum((sucs_peaks-mean(sucs_peaks)).^2);
        
        if R2_exp<R2_linear
            % set respones before steady state to nans
            steady_state=find(yf<((yf(1)-f.c)/4+f.c),1); % define steady state at point where 3/4 of the dacay has taken place
            temp = find(resp_idx,steady_state);
            steady_state = temp(end); % needs to be adjusted to mean steady state based on all trials instead of all sucesses
            if isempty(steady_state)
                steady_state = 1;
            end
        else
            steady_state=1;
        end
        
        resp_prop = sum(resp_idx(steady_state:end))/length(resp_idx(steady_state:end));
        temp = peaks(steady_state:end);
        avg_amp = mean(temp(find(resp_idx(steady_state:end))));
        
    catch
        steady_state = nan;
        resp_prop = sum(resp_idx)/length(resp_idx);
        avg_amp = mean(peaks(find(resp_idx)));
    end
    
else
    resp_idx = zeros(size(bs_traces,2),1);
    steady_state = nan;
    resp_prop = 0;
    avg_amp = 0;
end

end

function [Rs_cell, Rm_cell, Cm_cell]=calc_cell_acess_prop(header,bs_traces,sr,srF, base_end)
if header.ephys.ephys.stimOnArray==1
    
    testpulse_start =header.ephys.ephys.pulseParameters{1, 1}.squarePulseTrainDelay*sr;
    
    testpulse_end       =   testpulse_start+1000;
    
    
    test_window_ind=[(testpulse_start-1) : testpulse_end]';
    baseline_window_ind=[(1*srF+1) : base_end*srF]';
    
    numPointsTest=size(test_window_ind,1);
    numTraces=size(bs_traces,2);
    ydata=bs_traces;
    cellparam_ydata=zeros(numPointsTest, numTraces);
    try
        cellparam_ydata=ydata(test_window_ind,:);
    catch
        cellparam_ydata=ones(1000,numTraces)*NaN;
    end
    magdY=abs(diff(cellparam_ydata(:,1)));
    divs3=find(magdY>0.5.*max(magdY));
    
    baseline_mean=median(ydata(baseline_window_ind,:));
    
    try
        cellparam_ydata=cellparam_ydata(divs3(1)+round(srF./2):divs3(end)-round(srF./2),:);
    catch
        cellparam_ydata=cellparam_ydata;
    end
    Slope1=cellparam_ydata(1,:)-cellparam_ydata(2,:);
    Slope2=cellparam_ydata(2,:)-cellparam_ydata(3,:);
    dSlope=Slope1-Slope2;
    SlopeEst1=Slope1+dSlope;
    SlopeEst2=Slope1+2.*dSlope;
    
    Rs=0;
    Rm=0;
    Cm=0;
    amp=-5;
    pretestbase=baseline_mean;
    for j=1:numTraces
        midtestbase(j)=mean(cellparam_ydata(round(.3*end):end,j)');
        %get peak
        peak(j)=min(cellparam_ydata(:,j)');
        %improve estimate - first two points (sr=4 kHz) are missing!
        peak(j)=peak(j)+SlopeEst1(j).*1.3;
        %find tau
        r=find(cellparam_ydata(:,j) < (peak(j)-(midtestbase(j))).*(exp(-1))+midtestbase(j));
        tau(j) = length(r)./srF;   %time constant, in ms
        
        %%%%across holding potential and stepsper cell
        
    end
    
    Rs = abs(amp./(peak-pretestbase)); %in Gohms
    %Rm = abs(amp./(midtestbase-pretestbase)); %YP's code
    Rm = abs((amp-(midtestbase-pretestbase).*Rs)./(midtestbase-pretestbase)); %also in Gohms
    Cm=(Rs+Rm).*tau./(Rs.*Rm); %in pF
    Rs_cell=Rs*1000;
    Rm_cell=Rm*1000;
    Cm_cell=Cm;
    
else
    Rs_cell(size(bs_traces,2))=NaN;
    Rm_cell(size(bs_traces,2))=NaN;
    Cm_cell(size(bs_traces,2))=NaN;
end
end









