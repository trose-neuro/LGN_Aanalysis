function  [neg_failure, pos_failure]=minianalysis(list, idx, pathName, fc, show, user);
%SW181229
%Function that extracts minis by using the std threshold criterion

%list=      information of folders for each cell SW000XX
%idx=       vector with idx of which recording is mini recording
%pathName=  folder name of cell
%fc=        factor of how many stds the signal should be included
%show=      show plots or not (1 or 0)

base_start          =   1;
base_end            =   99;
%if user=0
pulse_start         =   100;
pulse_end           =   110;
%else
% blue_pulse_start         =   100;
% blue_pulse_end           =   110;

span = 1; %filter span
plotlength = 4000; %samples

%% filtering
filterephys = 1;        % filtering yes/no?
cutoff      = 1000       % Hz (use 500 Hz for mini event / amplitude detection. Chen & Regehr 2000)
order       = 4         % filter order ('pole')
type        = 'Bessel'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel)

if filterephys; 
    disp(['Filtering: ' num2str(order) ' pole ' type '-Filter w/ ' num2str(cutoff) ' Hz cutoff']);
end;
%% 

for i=1:length(idx);
    load([char(pathName) '/' list(idx(i)).name],'-mat');
    sr = header.ephys.ephys.sampleRate;%check sample rate
    srF = 1/(1000/sr);
    
    samples_per_sweep = header.ephys.ephys.traceLength*sr;
    
    ephystraces=data.ephys.trace_1;
    
    if filterephys
        ephystraces = lowpassfilt(ephystraces, order, cutoff, sr, type);
    end
    
    
    diodetraces=data.acquirer.trace_1;
    
    
    timebase=1:1/sr:samples_per_sweep/sr
    
    traces=reshape(ephystraces, samples_per_sweep, length(ephystraces)/samples_per_sweep); % why chunk it up in 1 sec bins hardcoded? Can the sweep length not be extracted from the header?
    diode=reshape(diodetraces, samples_per_sweep, length(ephystraces)/samples_per_sweep);
    
    bs=traces(base_start*srF:base_end*srF,:); %extract trace baseline periods
    diodebs=diode(base_start*srF:base_end*srF,:); %extract diodetrace baseline periods
    
    
    bs_std=std(bs);
    bs_traces=bsxfun(@minus, traces, mean(bs)); %mean-baseline-subtract each trial individually
    bs_diode=bsxfun(@minus, diode, mean(diodebs)); %mean-baseline-subtract each trial individually
    
    neg_peak=min(bs_traces(pulse_start*srF:pulse_end*srF+40*srF,:)); % why the 40 sample offset?
    pos_peak=max(bs_traces(pulse_start*srF:pulse_end*srF+40*srF,:));
    
    neg_fail=neg_peak<fc*bs_std*(-1); % threshold is set by z-score. Hmm... not sure
    pos_fail=pos_peak>fc*bs_std;
    
    neg_idx=find(neg_fail==1);
    pos_idx=find(pos_fail==1);
    
    neg_m=zeros(1,size(traces,2));
    neg_m(neg_idx)=neg_peak(neg_idx);
    pos_m=zeros(1,size(traces,2));
    pos_m(pos_idx)=pos_peak(pos_idx);
    
    try
        neg_fail(neg_idx)=neg_peak(neg_idx);
        pos_fail(pos_idx)=pos_peak(pos_idx);
    catch
        neg_fail=zeros(length(pos_peak));
        pos_fail=zeros(length(pos_peak));
    end
    neg_failure(:,i)=neg_m;
    pos_failure(:,i)=pos_m;
    
    %PLOT
    if show==1
        figure;
        for k=1:size(traces,2)
            subplot(size(bs_traces,2)/10,10,k);
            plot(bs_traces(1:plotlength,k), 'k'); hold on
            ylim([round(min(neg_peak)/10)*10 round(max(pos_peak)/10)*10]);
            hline(neg_peak(k));
            yyaxis right
            plot(smooth(bs_diode(1:plotlength,k),span)); ylim([-1 5]);
            yyaxis left
        end
        
        figure;
        for k=1:size(traces,2)
            ppp=plot(bs_traces(1:plotlength,k), '-k'); hold on
        end
        yyaxis right
        plot(smooth(bs_diode(1:plotlength,k),span)); %ylim([-1 5]);
        yyaxis left
        
        figure;
        for k=2
            ppp=plot(diff(bs_traces(1:plotlength,k)), '-k'); hold on
        end
       
        
        figure;
        plot(neg_peak,'*');
        hold on
        plot(pos_peak,'*');
        hold on;
        plot(repmat(fc*(mean(bs_std))*(-1),length(pos_peak)));
        hold on;
        plot(repmat(fc*(mean(bs_std)),length(pos_peak)));
    end
end
end
