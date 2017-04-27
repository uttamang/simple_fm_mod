% Hier fehlen noch
% Anpassung an f0. Aus gangfrequenz stimmt noch nicht.
% interpolation
clear all;
close all;
% gegebene Werte
f0 = 12*10^3;    % Grundfrequenz
del_f0 =4e3;     % 4kHz
f0min=f0-del_f0; % 8kHz
f0max=f0+del_f0; % 16kHz
fa = 48*10^3;    % Abtastfrequenz
f_tab = 10;      % Cos-Tabellen schritte
A = 32000;           % Amplitude
k= f0/32000;         %Frequenzhub
Ta = 1/fa;       % Abtastzeit
fx = 50; % Frequen des Eingangssignals
Tx = 1/fx; %Periodendauer des Eingangssignals
p=(2*pi/4800);   % cos-Tabelle-Schritte
n=1;
outp=16; % Ausgangsinformation in einer Abtastzeit
total_y_point=n*outp; % Ausgangspunkte in einer Abtastzeit
t_total = Tx/Ta; %  Zeit für eine Eingangssignalperiode
% Matrizen erstellen
cos_ta =  zeros (4800,1);        % cos-Tabelle 4800 Zeilen, 1 Spalten
x = zeros (2,t_total);             % Eingangssignal erst mal als Matrix mit 9600 Spalten
y = zeros (4,t_total*total_y_point); % Ausgangssignal als Matrix 9600*(Anzahl des Ausgangs für ein abgetastetes Signal)
%input_ping_buffer = zeros (32,1);
%input_pong_buffer = zeros (32,1);    

%output_ping_buffer = zeros (32,1);
%output_pong_buffer = zeros (32,1);

% cos-Tabelle erstellen
for ww =1:1:4800;
    cos_ta(ww,1) = cos(p*ww);    % berechnung cos Wert
    %cos_ta(w,2)= p*w;           % Winkel zum entsprechenden cos Wert eintragen
end;

% Eingangssignal erzeugen

for ee = 1:1:t_total;  
    x(1,ee) = A*sin(2*pi*ee/t_total);       % Eingangssignal erst mal ein sin mit 50 Hz
    % x(1,ee)=round(ee/480);
end;


del_ph_max=2*pi*Ta*(n*f0+k*max(x(1,:)));    % max-Phasendifferenz
del_ph_min=2*pi*Ta*(n*f0+k*min(x(1,:)));    % min-Phasendifferenz
for nn=1:1:t_total; % Durchlauf aller Abtaszeiten
    del_ph=2*pi*Ta*(n*f0+k*x(1,nn));    % Phasendifferenz berechnen
    ph_inc=(del_ph-(del_ph_min))/(del_ph_max-del_ph_min)*4800/outp*4+4800/outp; % entsprechende Phasenincrement berechnen
    ph_inc_ganz=round(ph_inc);          % Phasen increment aufrunden. Eigentlich muss hier durch eine Interpolartion ersetzen.
    %ph_inc_ganz_v(nn,1)=ph_inc_ganz;    % Phasenincrement abspeichern => Plot
    for ii=1:1:total_y_point;
        tt=(nn-1)*total_y_point+ii;     % Zeiger
        pp=ph_inc_ganz*tt;              % Zeiger Phasenincrement
        run=mod(pp,4800);               % Begrenzung des Zeigers
        y(1,tt)=cos_ta(run+1,1);        % Ausgangssignal aus der Tabelle holen
        y(2,tt)=del_ph;                 % Phasendifferenz ,nur für Plot
        y(3,tt)=ph_inc;                 % Phasendefferenz ,nur für Plot
    end;
end;
% Plots
ly=1:1:length(y);
yx=ly/outp;
f0_min=n*f0+k*min(x(1,:));
f0_max=n*f0+k*max(x(1,:));
subplot(2,1,1);
plot(x(1,1:t_total))
title('Eingangssignal x(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
legend(['max : ' num2str(max(x(1,:))),'  min : ' num2str(min(x(1,:)))]);
xlim([0 t_total]);
subplot(2,1,2);
plot(y(1,:))
title('Ausgangssignal y(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
xlim([0 length(y)]);

figure
subplot(4,1,1);
plot(y(1,:))
title('Ausgangssignal y(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
xlim([0 length(y)/4]);

subplot(4,1,2);
plot(y(1,:))
title('Ausgangssignal y(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
xlim([length(y)/4 length(y)/2]);

subplot(4,1,3);
plot(y(1,:))
title('Ausgangssignal y(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
xlim([length(y)/2 length(y)*3/4]);

subplot(4,1,4);
plot(y(1,:))
title('Ausgangssignal y(n)');
xlabel('Abtastpunkt n oder Zein in s/Ta');
xlim([length(y)*3/4 length(y)]);

figure
subplot(2,1,1);
plot(y(2,:))
title('delta-Phi');
xlabel('Abtastpunkt n oder Zein in s/Ta');
legend(['max : ' num2str(max(y(2,:))),'  min : ' num2str(min(y(2,:)))]);
xlim([0 length(y)]);

subplot(2,1,2);
plot(y(3,:))
xlim([0 length(y)]);
legend(['max : ' num2str(max(y(3,:))),'  min : ' num2str(min(y(3,:)))]);
title('Phasenincrement 4800 gespeicherten Tabelle');
xlabel('Abtastpunkt n oder Zein in s/Ta');
return;


for r = 1:1:300  ;              
    
    % Wechsel zwischen Ping- und Pong-Buffer
    if mod(r,2);                % odd
        n=1;
    else                        %even
        n=2;
    end;
        
    % Unterscheidung ob ping oder pong in Blöcken
    if n == 1;                  % wenn n == 1 dann ping verwenden
        for kk=1:1:32;           % Teil des Eingangssignals in "Buffer" schreiben dabei Spalten zu Zeilen
            input_ping_buffer(kk,1) = x(1,kk*r);
            del_Phi=2*pi*Ta+2*pi*1000*Ta*x(1,kk*r); % Phasendifferenz
        end;
    else                        % wenn n == 2 dann pong verwenden
        for kk=1:1:32;          
            input_pong_buffer(kk,1) = x(1,kk*r);
        end;
    end;
    
    % Berechnung für einen Input_Buffer durchlauf
    
figure(1);    
plot(cos_ta(:,1))
figure(2);   
plot(cos_ta(:,2))
figure(3);   
plot(x(1,:))

return;
% output_buffer mit cos-Wert aus der cos-Tabelle
    for z = 1:1:32;
            output_ping_buffer(1,z) = A*cos_ta(del_Phi);
    end;
    
    
    
    
    % Ausgangssignal von FM-mod
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
