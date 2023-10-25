% Confronto tra i vari metodi di discretizzazione 
%
s = tf('s');
gs = 25 / (s*(s+1) * (s+10));		  	% Funzione di trasferimento gs
k = 1; tau1 = 0.806; tau2 = 0.117;      % Parametri per rete correttrice
grs = k * ((1+tau1*s) / (1+tau2*s));    % Rete correttrice anticipatrice 




h=figure(1); clf;
hold on;
bode(grs,'b')                                                                              %graficazione di bode per tempo continuo

[mag_Gs,ph_Gs,wout]=bode(grs);                                                             %ricavati i valori di bode
grs_inf=k*tau1/tau2;                                                                       %valore all'infinito
children = get(h, 'Children');
magChild = children(3);                                                                    % Pick a handle to axes of magnitude in bode diagram
axes(magChild)
plot([10^-2 10^3],20*log10(grs_inf)*[1 1],'r')                                             %plottata la linea parallela al valore infinito
Ind=find(abs(mag_Gs(1,1,:)-grs_inf/sqrt(2))==min(abs(mag_Gs(1,1,:)-grs_inf/sqrt(2))));     %trovato indice al quale il valore di grs è inferiore di 3 db al valore all'infinito
plot(wout(Ind),20*log10(mag_Gs(1,1,Ind))*[1 1],'.m','MarkerSize',8)                        %plottato il punto corrispondente




%
for T=[0.01]	% Periodi di campionamento
    z=tf('z',T);
    %
    figure(2); hold off; clf
    Tfin=6;
    gsr=grs*gs;                         % Il nostro sistema finale, non ancora retroazionato
    gsr0=feedback(gsr,1);		        % Sistema retroazionato gs0=gs/(1+gs)
    [y,t]=step(gsr0,Tfin);		        % Risposta al gradino del sistema retroazionato per Tfin secondi
    plot(t,y,'b'); hold on		        % Graficazione con tratto di colore rosso di come reagisce il sistema tempo continuo al gradino unitario
    %
    gz=c2d(gs,T,'zoh');			        % Si discretizza (gz) il sistema gs
    %
    % 1) Differenze all'indietro
    %
    a=k*(T+tau1); b=-k*tau1; c=T+tau2; d=-tau2;
    grz1=(a*z+b)/(c*z+d)
    % 
    figure(1);
    bode(grz1,'r')                      %graficazione bode
    figure(2);
    %
    gz1=feedback(grz1*gz,1);		    % Sistama discreto retroazionato: gz0=gz/(1+gz)
    [y,t]=step(gz1,Tfin);		        % Risposta al gradino
    stairs(t,y,'r'); hold on		    % Graficazione gradino
    %
    % 2) Differenze all'avanti
    %
    a=k*tau1; b=k*(T-tau1); c=tau2; d=T-tau2;
    grz2=(a*z+b)/(c*z+d)
    %
    figure(1);
    bode(grz2,'g')                      %graficazione bode
    figure(2);
    %
    gz2=feedback(grz2*gz,1);		    % Sistama discreto retroazionato: gz0=gz/(1+gz)
    [y,t]=step(gz2,Tfin);		        % Risposta al gradico
    stairs(t,y,'g'); hold on		    % Graficazione gradino
    % Se graficato risulta instabile
    %
    % 3) Bilineare
    a=k*(T+ 2*tau1); b=k*(T-2*tau1); c=T+2*tau2; d=T-2*tau2;
    grz3=(a*z+b)/(c*z+d);
    %
    figure(1);
    bode(grz3,'m')                      % Graficazione Bode
    
    figure(2);
    %
    gz3=feedback(grz3*gz,1);		    % Sistama discreto retroazionato: gz0=gz/(1+gz)
    [y,t]=step(gz3,Tfin);		        % Risposta al gradico
    stairs(t,y,'m'); hold on		    % Graficazione gradino
    %
    % 4) Bilineare con precompensazione

    s=wout(Ind)/tan(wout(Ind)*T/2)*(1-z^-1)/(1+z^-1);   %fatta alla wout(Ind)
    
    grz4 = k * ((1+tau1*s) / (1+tau2*s));               % Rete correttrice anticipatrice 
    
    figure(1);
    bode(grz4,'y')                                      % Graficazione Bode
    figure(2);
    
    
    gz4 = feedback(grz4*gz,1);                          
    [y,t]=step(gz4,Tfin);
    stairs(t,y,'y'); hold on                            % Graficazione gradino

   % 7) corrispondenza poli-zeri
    alpha=exp(-T/tau1); beta=exp(-T/tau2)
    h=k*(1-beta)/(1-alpha);
    a=h; b=-h*alpha; c=1; d=-beta;
    grz7=(a*z+b)/(c*z+d)
    %
    figure(1);
    bode(grz7,'k')                                      % Graficazione Bode
    title('Gs (b), indietro (r), avanti (g), bilineare (m), bil comp (y), corr. pz (k)')
    %title('Gs (b), indietro (r), avanti (g), bilineare (m), bil comp (y)')
    figure(2);
    %
    gz7=feedback(grz7*gz,1);		                    % Sistama discreto retroazionato: gz0=gz/(1+gz)
    [y,t]=step(gz7,Tfin);		                        % Risposta al gradino
    stairs(t,y,'k'); hold on		                    % Graficazione 
    axis([0 Tfin 0 1.5])
    xlabel('Tempo (s)')
    ylabel('Im[G(j*w)]')
    title('T cont (b), indietro (r), avanti (g), bilineare (m), bil comp (y), corr. pz (k)')
    %title('T cont (b), indietro (r), avanti (g), bilineare (m), bil comp (y)')

    pause
end

%
return%
