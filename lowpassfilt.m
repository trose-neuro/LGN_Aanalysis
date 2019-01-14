function filtered = lowpassfilt(signal, order, fc, fs, type)

%TR2019 switchboard to standard ML 'analog' filters.

if nargin < 5
    type = 'Bessel';
end

if type == 'Butter'
    [b,a] = butter(order,fc/(fs/2), 'low');
elseif type == 'Bessel'
    [b,a] = besself(order,fc * 2 * pi(), 'low');
    [b,a] = bilinear(b,a,fs); %digital to analog conversion (besself is analog only)
end

filtered = filtfilt(b,a,signal); % bidirectional filtering w/o phase distortion

  