%% constants
close all 
clear all
pi = 3.1415926;
debug = 1; % debug mode, show several plotted matrices

%% frequencies
fs = 48000; % sampling frequency
fc = 12000; % center frequency

Ts = 1/fs;
f_sig = 680; % signal frequency in Hz

A = 1; % signal amplitude
signal_time_step = Ts;
t = 0.1; % signal sample length in sigure;
time_array = (0:signal_time_step:t-signal_time_step);

lut_table_size = 4800; % length for the cosine lookup table
buffer_length = 64; % length of ping and pong buffers
nk=3;  % 8 waere besser
k = fc/nk; % frequency deviation 
f_out_max = fc+k;
f_out_min = fc-k;
lut_fak=300/(1-1/nk);

%% generate signal
signal = A*cos(2*pi*f_sig*time_array-pi/2);
% signal = A*cos(2*pi.*rand()*500*time_array-pi/2);

%% generate lut
cos_lut = cos(2*pi*(1:1:lut_table_size)/lut_table_size);

if debug == 1;
figure;
plot(1:1:lut_table_size,cos_lut)
title('Cosine lookup table');
end;

%% samples in one buffer
n_samples = 16;
% f_output =?
signal_periods_in_buffer = buffer_length/n_samples;

if debug == 1;
figure;
plot(time_array, signal);
title('Plotted input signal');
end;
%% read input signal to buffer

inputbuffer = zeros (2,buffer_length);
outputbuffer = zeros (2,buffer_length*n_samples);
buffer = 1;
pointer_old=0;
 for r = 1:1:20;
    
     for i = 1:1:buffer_length;
         n = (r-1)*buffer_length+i; % sample time 
         inputbuffer(buffer,i) = signal(n);  % ping -> pong
         for j = 1:1:n_samples;
            kk =  (i-1)*n_samples+j;
            del_phi = 2*pi*(1/fs)*fc+2*pi*(1/fs)*k*inputbuffer(buffer,i);
            pointer_new=del_phi/(2*pi)*4*lut_fak;
            pointer_min = lut_fak*(1-1/nk);
            pointer_corr = pointer_new+pointer_old;
            pointer_corr = round(pointer_corr);
            pointer_corr = mod(pointer_corr, lut_table_size);
            pointer_old=pointer_corr;
            outputbuffer(buffer,kk) = cos_lut(pointer_corr+1);
        end;
     end;
    
     buffer = not(buffer-1)+1; % change buffer from ping to pong or the other way
   
 end;

figure;
subplot(3,1,1)
plot(0:1:buffer_length-1,inputbuffer(1,:));
hold on 
plot(0:1:buffer_length-1,inputbuffer(2,:));
xlim([0 buffer_length-1]);
legend('ping-in','pong-in');
title('last input buffer contents')
subplot(3,1,2)
plot(outputbuffer(1,:));
xlim([0 length(outputbuffer(1,:))]);
legend('ping-out');
title(['last output buffer contents f-out = [' ,num2str(f_out_min) ':' num2str(f_out_max) ']Hz'])
subplot(3,1,3)
plot(outputbuffer(2,:),'r');
xlim([0 length(outputbuffer(2,:))]);
legend('pong-out');
title(['last output buffer contents f-out = [' ,num2str(f_out_min) ':' num2str(f_out_max) ']Hz'])
