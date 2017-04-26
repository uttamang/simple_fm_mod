% asdj
% constants

f0 = 12*10^3;    % centerfrequency
fa = 48*10^3;    % samplingfrequency
f_tab = 10;      % cos table step size
A = 1;           % amplitude
k= 1000;         % frequency deviation
 
Ta = 1/fa;       % sampling time
p=(2*pi/4800);   % cos-sampling steps

% create matrices
cos_ta = zeros (4800,2);       
x = zeros (1,9600);             % inputsignal

input_ping_buffer = zeros (32,1);
input_pong_buffer = zeros (32,1);    

output_ping_buffer = zeros (32,1);
output_pong_buffer = zeros (32,1);

y=zeros (1,9600);

% initialize cos table
q=0;
for w =1:1:4800;
    cos_ta(w,1) = cos(p*q);    % calculate cos values
    cos_ta(w,2)=p*q;           % phase angle values
    q=q+1;
end;

% initialize input signal

for e = 1:1:9600;
    x(1,e) = sin(2*pi*e/9600);     % inputsignal: sin 50Hz
end;

% 300 mal 32er teile in buffer schreiben -> nur noetig wenn wir begrenzt grosses Eingangssignal haben

for r = 1:1:300  ;              
    
    % change ping and pong buffer
    if mod(r,2);                % odd
        n=1;
    else                        % even
        n=2;
    end;
        
    if n == 1;                  % use ping
        for k=1:1:32;           % write a section of the input signal
            input_ping_buffer(k,1) = x(1,k*r);
        end;
    else                        % use pong
        for k=1:1:32;          
            input_pong_buffer(k,1) = x(1,k*r);
        end;
    end;
    
    % calculation after one input buffer is filled
    
    
	% output buffer
    for z = 1:1:32;
            output_ping_buffer(1,z) = A*cos_ta(m);
    end;
    
    
    
    
    % outputsignal of the fm modulation
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
