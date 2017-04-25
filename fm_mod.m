
% gegebene Werte

f0 = 12*10^3;    % Grundfrequenz
fa = 48*10^3;    % Abtastfrequenz
f_tab = 10;      % Cos-Tabellen schritte
A = 1;           % Amplitude
k= 1000;         %Frequenzhub
 
Ta = 1/fa;       % Abtastzeit
p=(2*pi/4800);   % cos-Tabelle-Schritte

% Matrizen erstellen
cos_ta = zeros (4800,2);        %4800 Zeilen, 2 Spalten
x = zeros (1,9600);             % Eingangssignal erst mal als Matrix mit 9600 Spalten

input_ping_buffer = zeros (32,1);
input_pong_buffer = zeros (32,1);    

output_ping_buffer = zeros (32,1);
output_pong_buffer = zeros (32,1);

y=zeros (1,9600);

% cos-Tabelle erstellen
q=0;
for w =1:1:4800;
    cos_ta(w,1) = cos(p*q);    % berechnung cos Wert
    cos_ta(w,2)=p*q;           % Winkel zum entsprechenden cos Wert eintragen
    q=q+1;
end;

% Eingangssignal erzeugen

for e = 1:1:9600;
    x(1,e) = sin(2*pi*e/9600);     % Eingangssignal erst mal ein sin mit 50 Hz
end;

% 300 mal 32er teile in buffer schreiben -> nur noetig wenn wir begrenzt grosses Eingangssignal haben

for r = 1:1:300  ;              
    
    % Wechsel zwischen Ping- und Pong-Buffer
    if mod(r,2);                % odd
        n=1;
    else                        %even
        n=2;
    end;
        
    % Unterscheidung ob ping oder pong in Bl�cken
    if n == 1;                  % wenn n == 1 dann ping verwenden
        for k=1:1:32;           % Teil des Eingangssignals in "Buffer" schreiben dabei Spalten zu Zeilen
            input_ping_buffer(k,1) = x(1,k*r);
        end;
    else                        % wenn n == 2 dann pong verwenden
        for k=1:1:32;          
            input_pong_buffer(k,1) = x(1,k*r);
        end;
    end;
    
    % Berechnung f�r einen Input_Buffer durchlauf
    
    
    
    
    
	 % output_buffer mit cos-Wert aus der cos-Tabelle
    for z = 1:1:32;
            output_ping_buffer(1,z) = A*cos_ta(m);
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
