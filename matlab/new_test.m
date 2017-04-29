%% constants
pi = 3.1415926;
debug = 0; % debug mode, show several plotted matrices

fs = 48000; % sampling frequency
fc = 12000; % center frequency

Ts = 1/fs;

A = 1; % signal amplitude
signal_time_step = Ts;
t = 4; % signal sample length in sigure;
time_array = (0:signal_time_step:t-signal_time_step);
f_sig = 7000; % signal frequency in Hz
lut_table_size = 4800; % length for the cosine lookup table
buffer_length = 32; % length of ping and pong buffers

%% generate signal
signal = A*cos(2*pi*f_sig*time_array-pi/2);

%% generate lut
cos_lut = cos(2*pi*(1:1:lut_table_size)/lut_table_size);

if debug == 1;
figure;
plot(1:1:lut_table_size,cos_lut)
title('Cosine lookup table');
end;

%% samples in one buffer
n_samples = 1/(Ts*f_sig);
signal_periods_in_buffer = buffer_length/n_samples;

if debug == 1;
figure;
plot(time_array, signal);
title('Plotted input signal');
end;

%% read input signal to buffer

buffer = [1:1:buffer_length;1:1:buffer_length];
n = 1;
for r = 1:1:20;
    n = not(n-1)+1;
    for i = 1:1:buffer_length;
        buffer(n,i) = signal(r*32+i);
    end;
end;

figure;
subplot(2,1,1)
plot(0:1:31,buffer(1,1:1:32));
title('last input ping buffer contents')
subplot(2,1,2)
plot(0:1:31,buffer(2,1:1:32));
title('last input pong buffer contents')
