function filtered = butter_lowpassfilt(signal, order, fc, fs)

[b,a] = butter(order,fc/(fs/2));

filtered = filter(b,a,signal);
