filename = '~/data/SW0003/SW0003AAAA0051.xsg'

load(filename, '-mat')
ephystraces=data.ephys.trace_1;

% ephystraces = ephystraces - mean(ephystraces(1:1000)) ;



%% TR2019: filtering
% filterephys = 1;        % filtering yes/no?
cutoff      = 5000;      % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Bessel'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel. Use Bessel at > 4 order to prevent ripples)
p           = 1 ; %detrending polynomial


sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
samples_per_sweep = header.ephys.ephys.traceLength*sr;


offset = samples_per_sweep/2;
timebase_all=offset/sr:1/sr:size(ephystraces,1)/sr; %TR2019: timebase
timebase=offset/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase

% ephystraces = lowpassfilt(ephystraces, order, cutoff, sr, type);


traces=reshape(ephystraces, samples_per_sweep, length(ephystraces)/samples_per_sweep)';

traces = traces(1,offset:end);

% save data to runresults
filename = '~/data/test-data.mat';
infername = '~/data/test-data-0001.mat';
save(filename,'traces')


%% extraction parameters
% build a parameter struct that will point to this file


params.burn_in_sweeps = 50;
params.num_sweeps = 100;
params.par = 0;
params.traces_filename = filename;
params.event_sign = -1;
params.dt = 1/sr;
params.exclusion_bound = 10/sr;
params.tau1_prop_std = 2/sr;
params.tau2_prop_std = 20/sr;
params.a_min = 5;
params.a_max = 500;
% fill rest of struct with defaults
params = get_params(params);

% check that you're not writing over a previous results file
if exist(params.full_save_string, 'file') == 2
    disp('****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****')
    disp(params.full_save_string)
    return
end

% run sampler
tic
run_posterior_sampler(params);
toc

load(filename, '-mat')
load(infername, '-mat')

try
    clear trialtimes amplitudes
end

close all

for i = params.burn_in_sweeps:size(results.trials.base,2);
    
    
    trialtimes{i} = results.trials.times(sum( results.trials.num_events(1:i-1))+1:...
        sum( results.trials.num_events(1:i-1))+ results.trials.num_events(i));
    
    amplitudes{i} = results.trials.amp(sum( results.trials.num_events(1:i-1))+1:...
        sum( results.trials.num_events(1:i-1))+ results.trials.num_events(i));
    
    tau_1{i} = results.trials.tau1(sum( results.trials.num_events(1:i-1))+1:...
        sum( results.trials.num_events(1:i-1))+ results.trials.num_events(i));
    
    tau_2{i} = results.trials.tau2(sum( results.trials.num_events(1:i-1))+1:...
        sum( results.trials.num_events(1:i-1))+ results.trials.num_events(i));
    
    
end

figure(43625)

subplot(3,1,1), plot(traces);hold all


subplot(3,1,2),
[p k ]= hist([trialtimes{:}],samples_per_sweep/sr*1000);ylim('auto'); hold all
plot(k,p);
[peakamppeak_cnt peakloc] = findpeaks(p,k,'MinPeakDistance',100, 'MinPeakHeight', std(p));
plot(peakloc, peakamppeak_cnt, 'ok');

subplot(3,1,3),
[pa ka ]= hist([amplitudes{:}],samples_per_sweep/sr*1000);ylim('auto'); hold all
plot(ka,pa);

% [peakamppeak_cnt peakloc] = findpeaks(p,k,'MinPeakDistance',100, 'MinPeakHeight', std(p));
% plot(peakloc, peakamppeak_cnt, 'ok');


subplot(3,1,1), vline(peakloc);

% plot(peakloc, [amplitudes(peakloc)*-1, 'ok');


% subplot(3,1,3),
% hist([amplitudes{1000:2000}],1000);%ylim([0 1000])


