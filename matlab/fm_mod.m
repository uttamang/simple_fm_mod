% Hier fehlen noch
% Anpassung an f0. Aus gangfrequenz stimmt noch nicht.
% interpolation

clear all;
close all;
% constants
f0 = 12*10^3;    % center frequency
del_f0 =4e3;     % 
f0min=f0-del_f0; % min. modulation frequency =  8kHz
f0max=f0+del_f0; % max. modulation frequency = 16kHz
fa = 48*10^3;    % sampling frequency
f_tab = 10;      % cos lut step size
A = 32000;       % amplitude
k= f0/32000;     % frequency deviation
Ta = 1/fa;       % sampling time
fx = 50;         % input signal frequency
Tx = 1/fx;       % input signal period
p=(2*pi/4800);   % cos lut steps
n=1;
outp=16;         % Ausgangsinformation in einer Abtastzeit
total_y_point=n*outp; % Ausgangspunkte in einer Abtastzeit
t_total = Tx/Ta; %  Zeit f�r eine Eingangssignalperiode
% Matrizen erstellen
cos_ta =  zeros (4800,1);            % fill cos lut table with zeroes
x = zeros (2,t_total);               % fill input signal matrix with zeroes
y = zeros (4,t_total*total_y_point); % fill output signal matrix with zeroes

%input_ping_buffer = zeros (32,1);
%input_pong_buffer = zeros (32,1);    

%output_ping_buffer = zeros (32,1);
%output_pong_buffer = zeros (32,1);

% create cos table
for ww =1:1:4800;
    cos_ta(ww,1) = cos(p*ww);    % calculate cos values
    %cos_ta(w,2)= p*w;           % put phase to cos values
end;

% create input signal

for ee = 1:1:t_total;  
    x(1,ee) = A*sin(2*pi*ee/t_total);       % 
    % x(1,ee)=round(ee/480);
end;


del_ph_max=2*pi*Ta*(n*f0+k*max(x(1,:)));    % max-phase offset
del_ph_min=2*pi*Ta*(n*f0+k*min(x(1,:)));    % min-phase offset
for nn=1:1:t_total; % run through all sampling times
    del_ph=2*pi*Ta*(n*f0+k*x(1,nn));     % calculate phase offset
    ph_inc=(del_ph-(del_ph_min))/(del_ph_max-del_ph_min)*4800/outp*4+4800/outp; % calculate phase increment
    ph_inc_ganz=round(ph_inc);           % round phase increment, replace with linear interpolation
    %ph_inc_ganz_v(nn,1)=ph_inc_ganz;    % save phase increment
    for ii=1:1:total_y_point;
        tt=(nn-1)*total_y_point+ii;     % pointer
        pp=ph_inc_ganz*tt;              % pointer phase increment
        run=mod(pp,4800);               % wrap pointer
        y(1,tt)=cos_ta(run+1,1);        % get outputsignal from table
        y(2,tt)=del_ph;                 % phaseoffset for plot
        y(3,tt)=ph_inc;                 % phaseoffset for plot
    end;
end;


% phlots
ly=1:1:length(y);
yx=ly/outp;
f0_min=n*f0+k*min(x(1,:));
f0_max=n*f0+k*max(x(1,:));
subplot(2,1,1);
plot(x(1,1:t_total))
title('input signal x(n)');
xlabel('sample point n or z(in) in s/Ta');
legend(['max : ' num2str(max(x(1,:))),'  min : ' num2str(min(x(1,:)))]);
xlim([0 t_total]);
subplot(2,1,2);
plot(y(1,:))
title('output signal y(n)');
xlabel('sample n oder z(in) in s/Ta');
xlim([0 length(y)]);

figure
subplot(4,1,1);
plot(y(1,:))
title('output signal y(n)');
xlabel('sample n oder Z(in) in s/Ta');
xlim([0 length(y)/4]);

subplot(4,1,2);
plot(y(1,:))
title('output signal y(n)');
xlabel('sample n oder Z(in) in s/Ta');
xlim([length(y)/4 length(y)/2]);

subplot(4,1,3);
plot(y(1,:))
title('output signal y(n)');
xlabel('sample n oder Z(in) in s/Ta');
xlim([length(y)/2 length(y)*3/4]);

subplot(4,1,4);
plot(y(1,:))
title('output signal y(n)');
xlabel('sample n oder Z(in) in s/Ta');
xlim([length(y)*3/4 length(y)]);

figure
subplot(2,1,1);
plot(y(2,:))
title('delta-phi');
xlabel('sample n oder Z(in) in s/Ta');
legend(['max : ' num2str(max(y(2,:))),'  min : ' num2str(min(y(2,:)))]);
xlim([0 length(y)]);

subplot(2,1,2);
plot(y(3,:))
xlim([0 length(y)]);
legend(['max : ' num2str(max(y(3,:))),'  min : ' num2str(min(y(3,:)))]);
title('phasen increment');
xlabel('sample n oder Z(in) in s/Ta');
return;


for r = 1:1:300  ;              
    
    % change buffer
    if mod(r,2); % odd
        n=0;
    else         % even
        n=1;
    end;
        
    % decide between ping or pong
    if n == 0;                  % if 0 then use ping
        for kk=1:1:32;           % Teil des Eingangssignals in "Buffer" schreiben dabei Spalten zu Zeilen
            input_ping_buffer(kk,1) = x(1,kk*r);
            del_Phi=2*pi*Ta+2*pi*1000*Ta*x(1,kk*r); % Phasendifferenz
        end;
    else                        % if n == 1 use pong
        for kk=1:1:32;          
            input_pong_buffer(kk,1) = x(1,kk*r);
        end;
    end;
    
    % calculation through one buffer
    
figure(1);    
plot(cos_ta(:,1))
figure(2);   
plot(cos_ta(:,2))
figure(3);   
plot(x(1,:))

return;
% output_buffer with cos values from lut
    for z = 1:1:32;
            output_ping_buffer(1,z) = A*cos_ta(del_Phi);
    end;
    
    
    
    
    % outputsignal / fm signal
    if n ==1; 
        for u=1:1:32;
            y(1,u*r) = output_ping_buffer(u,1);
        end;      
    else

        for u=1:1:32;
            y(1,u*r) = output_pong_buffer(u,1);
        end;
    end;
    
end; 
