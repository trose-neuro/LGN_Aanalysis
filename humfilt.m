function filtered = humfilt(signal, order, fc, fs)
%LOWPASSFILT  Switchboard for ML bandpass filter (Butterworth)
%   filtered = LOWPASSFILT(signal, order, fc, fs, type) filters the input
%   signal at a notch freqeuency (fc, in Hz plusminus 1 Hz) of sampling rate (fs) using a
%   Butterworth (type = 'Butter') filter of 'order' order (aka pole).
%   see https://de.mathworks.com/help/signal/ug/remove-the-60-hz-hum-from-a-signal.html
%   CAVE! This is not recommended! Adaptive pre-digitization filtering
%   (humbug) is better (no ring etc). Even better: get rid of the hum!
% Tobias Rose 2019


d = designfilt('bandstopiir','FilterOrder',order, ...
               'HalfPowerFrequency1',fc-1,'HalfPowerFrequency2',fc+1, ...
               'DesignMethod','butter','SampleRate',fs);

filtered = filtfilt(d,signal);
