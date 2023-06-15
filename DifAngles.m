clear all;clc;

N=10000;
Alpha=0.002382;             
lambda=120;                 
kappa=2.645;                
bp=0.771; jp=148.16;         
alpha = 0.0025;             
rho = 0.76;                  
y0=1000;                     
bm=0.8;Bm=1.041231831;  
%%% kellogg 30

Emin=50;
Emax=1700;
Theta_min=25.9;
Theta_max=34.1;  
N=10000;


Alpha=0.002382;              % constant A
lambda=120;                  % absorption mean free path 120 g/cm2
kappa=2.645;                 % exponent (-)
bp=0.771; jp=148.16;         % correction factor (-); factor (GeV)
alpha = 0.0025;              % muon energy loss in GeV/g/cm2
rho = 0.76;                  % fraction of pion energy that is transferred to muon
y0=1000;                     % atmoshperic depth g/cm2
bm=0.8;Bm=1.041231831;       % correction factor (-); 

theta1=Theta_min*pi/180;theta2=Theta_max*pi/180;             % change from degrees to radians
for i=1:N
tm=theta1+(theta2-theta1).*rand();
tm1=2*tm/(pi);
theta_corr(i)=acos((1-tm1)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2=[theta1:0.01:theta2];                                % zenith angle lower and upper limits
theta3=cos(tm2).^2;                                      % angular distribution modelled as cosine squared
C4=max(theta3);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1=(theta1+(theta2-theta1)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    u=rand();                                            % select a random number from uniform U(0,1)
    f1=cos(Theta_m1)^2;                                  % calculate random variable X             
    f2=C4;                                            
    f3=f1/f2;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if u<=f3                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Theta(i)=Theta_m1;                          % accept zenith angle
        break
    end
end

Em=[Emin:0.1:Emax];                                                 % energy lower and upper limits
Ep=(Em+alpha*y0*(sec(Muon_Theta(i))-0.1))*(1/rho);
Pm=(0.1*cos(Muon_Theta(i)).*(1-(alpha*(y0*sec(Muon_Theta(i))-100)./(rho.*Ep)))).^(Bm./((rho.*Ep+100*alpha)*cos(Muon_Theta(i))));
f=Alpha.*Pm.*(Ep.^(-kappa))*lambda*bp*jp./(Ep.*cos(Muon_Theta(i))+bp*jp);
C6=max(f);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1=(Emin+(Emax-Emin)*rand());                           % select a random energy from uniform U(Emin,Emax)
    u=rand();                                                % select a random number from uniform U(0,1)
    Ep1=(Em1+alpha*y0*(sec(Muon_Theta(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1=(0.1*cos(Muon_Theta(i))*(1-(alpha*(y0*sec(Muon_Theta(i))-100)/(rho*Ep1))))^(Bm/((rho*Ep1+100*alpha)*cos(Muon_Theta(i))));   % calculate probability
    f4=Alpha*Pm1*(Ep1^(-kappa))*lambda*bp*jp/(Ep1*cos(Muon_Theta(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5=C6;                                            
    f6=f4/f5;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if u<=f6                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energy(i)=Em1;                                   % accept energy
        break
    end
end
end

Muon_table=[Muon_Theta;theta_corr;Muon_Energy]';

bin1=2*(numel(Muon_Energy)^(1/3));      % bin size for energies selected according to Rice rule
[C,x2]=hist(Muon_Energy,bin1);
C1=C/trapz(x2,C);                       % normalization of the histogram
bin8=2*(numel(theta_corr)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8,x8]=hist(theta_corr,bin8);
C9=C8/trapz(x8,C8);                     % normalzation of the histogram
bin10=2*(numel(Muon_Theta)^(1/3));      % bin size for angles selected according to Rice rule
[C10,x10]=hist(Muon_Theta,bin10);
C11=C10/trapz(x10,C10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kellogg 75

Eminkell=50;
Emaxkell=1700;
Theta_minkell=70.9;
Theta_maxkell=79.1;  

theta1kell=Theta_minkell*pi/180; theta2kell=Theta_maxkell*pi/180;             % change from degrees to radians
for i=1:N
tmkell=theta1kell+(theta2kell-theta1kell).*rand();
tm1kell=2*tmkell/(pi);
theta_corrkell(i)=acos((1-tm1kell)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2kell=[theta1kell:0.01:theta2kell];                                % zenith angle lower and upper limits
theta3kell=cos(tm2kell).^2;                                      % angular distribution modelled as cosine squared
C4kell=max(theta3kell);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1kell=(theta1kell+(theta2kell-theta1kell)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    ukell=rand();                                            % select a random number from uniform U(0,1)
    f1kell=cos(Theta_m1kell)^2;                                  % calculate random variable X             
    f2kell=C4kell;                                            
    f3kell=f1kell/f2kell;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ukell<=f3kell                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Thetakell(i)=Theta_m1kell;                          % accept zenith angle
        break
    end
end

Emkell=[Eminkell:0.1:Emaxkell];                                                 % energy lower and upper limits
Epkell=(Emkell+alpha*y0*(sec(Muon_Thetakell(i))-0.1))*(1/rho);
Pmkell=(0.1*cos(Muon_Thetakell(i)).*(1-(alpha*(y0*sec(Muon_Thetakell(i))-100)./(rho.*Epkell)))).^(Bm./((rho.*Epkell+100*alpha)*cos(Muon_Thetakell(i))));
fkell=Alpha.*Pmkell.*(Epkell.^(-kappa))*lambda*bp*jp./(Epkell.*cos(Muon_Thetakell(i))+bp*jp);
C6kell=max(fkell);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1kell=(Eminkell+(Emaxkell-Eminkell)*rand());                           % select a random energy from uniform U(Emin,Emax)
    ukell=rand();                                                % select a random number from uniform U(0,1)
    Ep1kell=(Em1kell+alpha*y0*(sec(Muon_Thetakell(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1kell=(0.1*cos(Muon_Thetakell(i))*(1-(alpha*(y0*sec(Muon_Thetakell(i))-100)/(rho*Ep1kell))))^(Bm/((rho*Ep1kell+100*alpha)*cos(Muon_Thetakell(i))));   % calculate probability
    f4kell=Alpha*Pm1kell*(Ep1kell^(-kappa))*lambda*bp*jp/(Ep1kell*cos(Muon_Thetakell(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5kell=C6kell;                                            
    f6kell=f4kell/f5kell;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ukell<=f6kell                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energykell(i)=Em1kell;                                   % accept energy
        break
    end
end
end

Muon_tablekell=[Muon_Thetakell;theta_corrkell;Muon_Energykell]';

bin1kell=2*(numel(Muon_Energykell)^(1/3));      % bin size for energies selected according to Rice rule
[Ckell,x2kell]=hist(Muon_Energykell,bin1kell);
C1kell=Ckell/trapz(x2kell,Ckell);                       % normalization of the histogram
bin8kell=2*(numel(theta_corrkell)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8kell,x8kell]=hist(theta_corrkell,bin8kell);
C9kell=C8kell/trapz(x8kell,C8kell);                     % normalzation of the histogram
bin10kell=2*(numel(Muon_Thetakell)^(1/3));      % bin size for angles selected according to Rice rule
[C10kell,x10kell]=hist(Muon_Thetakell,bin10kell);
C11kell=C10kell/trapz(x10kell,C10kell);  

%%%%%%%%%%%%%%%%%%%%%%%% Tsuji 30
Emints30=2;
Emaxts30=250;
Theta_mints30=26;
Theta_maxts30=34;  

theta1ts30=Theta_mints30*pi/180; theta2ts30=Theta_maxts30*pi/180;             % change from degrees to radians
for i=1:N
tmts30=theta1ts30+(theta2ts30-theta1ts30).*rand();
tm1ts30=2*tmts30/(pi);
theta_corrts30(i)=acos((1-tm1ts30)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2ts30=[theta1ts30:0.01:theta2ts30];                                % zenith angle lower and upper limits
theta3ts30=cos(tm2ts30).^2;                                      % angular distribution modelled as cosine squared
C4ts30=max(theta3ts30);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1ts30=(theta1ts30+(theta2ts30-theta1ts30)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    uts30=rand();                                            % select a random number from uniform U(0,1)
    f1ts30=cos(Theta_m1ts30)^2;                                  % calculate random variable X             
    f2ts30=C4ts30;                                            
    f3ts30=f1ts30/f2ts30;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uts30<=f3ts30                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Thetats30(i)=Theta_m1ts30;                          % accept zenith angle
        break
    end
end

Emts30=[Emints30:0.1:Emaxts30];                                                 % energy lower and upper limits
Epts30=(Emts30+alpha*y0*(sec(Muon_Thetats30(i))-0.1))*(1/rho);
Pmts30=(0.1*cos(Muon_Thetats30(i)).*(1-(alpha*(y0*sec(Muon_Thetats30(i))-100)./(rho.*Epts30)))).^(Bm./((rho.*Epts30+100*alpha)*cos(Muon_Thetats30(i))));
fts30=Alpha.*Pmts30.*(Epts30.^(-kappa))*lambda*bp*jp./(Epts30.*cos(Muon_Thetats30(i))+bp*jp);
C6ts30=max(fts30);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1ts30=(Emints30+(Emaxts30-Emints30)*rand());                           % select a random energy from uniform U(Emin,Emax)
    uts30=rand();                                                % select a random number from uniform U(0,1)
    Ep1ts30=(Em1ts30+alpha*y0*(sec(Muon_Thetats30(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1ts30=(0.1*cos(Muon_Thetats30(i))*(1-(alpha*(y0*sec(Muon_Thetats30(i))-100)/(rho*Ep1ts30))))^(Bm/((rho*Ep1ts30+100*alpha)*cos(Muon_Thetats30(i))));   % calculate probability
    f4ts30=Alpha*Pm1ts30*(Ep1ts30^(-kappa))*lambda*bp*jp/(Ep1ts30*cos(Muon_Thetats30(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5ts30=C6ts30;                                            
    f6ts30=f4ts30/f5ts30;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uts30<=f6ts30                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energyts30(i)=Em1ts30;                                   % accept energy
        break
    end
end
end

Muon_tablets30=[Muon_Thetats30;theta_corrts30;Muon_Energyts30]';

bin1ts30=2*(numel(Muon_Energyts30)^(1/3));      % bin size for energies selected according to Rice rule
[Cts30,x2ts30]=hist(Muon_Energyts30,bin1ts30);
C1ts30=Cts30/trapz(x2ts30,Cts30);                       % normalization of the histogram
bin8ts30=2*(numel(theta_corrts30)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8ts30,x8ts30]=hist(theta_corrts30,bin8ts30);
C9ts30=C8ts30/trapz(x8ts30,C8ts30);                     % normalzation of the histogram
bin10ts30=2*(numel(Muon_Thetats30)^(1/3));      % bin size for angles selected according to Rice rule
[C10ts30,x10ts30]=hist(Muon_Thetats30,bin10ts30);
C11ts30=C10ts30/trapz(x10ts30,C10ts30);  

%%%%%%%%%%%%%%% jokisch
Eminjo=1;
Emaxjo=1000;
Theta_minjo=68;
Theta_maxjo=82;  

theta1jo=Theta_minjo*pi/180; theta2jo=Theta_maxjo*pi/180;             % change from degrees to radians
for i=1:N
tmjo=theta1jo+(theta2jo-theta1jo).*rand();
tm1jo=2*tmjo/(pi);
theta_corrjo(i)=acos((1-tm1jo)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2jo=[theta1jo:0.01:theta2jo];                                % zenith angle lower and upper limits
theta3jo=cos(tm2jo).^2;                                      % angular distribution modelled as cosine squared
C4jo=max(theta3jo);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1jo=(theta1jo+(theta2jo-theta1jo)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    ujo=rand();                                            % select a random number from uniform U(0,1)
    f1jo=cos(Theta_m1jo)^2;                                  % calculate random variable X             
    f2jo=C4jo;                                            
    f3jo=f1jo/f2jo;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ujo<=f3jo                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Thetajo(i)=Theta_m1jo;                          % accept zenith angle
        break
    end
end

Emjo=[Eminjo:0.1:Emaxjo];                                                 % energy lower and upper limits
Epjo=(Emjo+alpha*y0*(sec(Muon_Thetajo(i))-0.1))*(1/rho);
Pmjo=(0.1*cos(Muon_Thetajo(i)).*(1-(alpha*(y0*sec(Muon_Thetajo(i))-100)./(rho.*Epjo)))).^(Bm./((rho.*Epjo+100*alpha)*cos(Muon_Thetajo(i))));
fjo=Alpha.*Pmjo.*(Epjo.^(-kappa))*lambda*bp*jp./(Epjo.*cos(Muon_Thetajo(i))+bp*jp);
C6jo=max(fjo);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1jo=(Eminjo+(Emaxjo-Eminjo)*rand());                           % select a random energy from uniform U(Emin,Emax)
    ujo=rand();                                                % select a random number from uniform U(0,1)
    Ep1jo=(Em1jo+alpha*y0*(sec(Muon_Thetajo(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1jo=(0.1*cos(Muon_Thetajo(i))*(1-(alpha*(y0*sec(Muon_Thetajo(i))-100)/(rho*Ep1jo))))^(Bm/((rho*Ep1jo+100*alpha)*cos(Muon_Thetajo(i))));   % calculate probability
    f4jo=Alpha*Pm1jo*(Ep1jo^(-kappa))*lambda*bp*jp/(Ep1jo*cos(Muon_Thetajo(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5jo=C6jo;                                            
    f6jo=f4jo/f5jo;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ujo<=f6jo                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energyjo(i)=Em1jo;                                   % accept energy
        break
    end
end
end

Muon_tablejo=[Muon_Thetajo;theta_corrjo;Muon_Energyjo]';

bin1jo=2*(numel(Muon_Energyjo)^(1/3));      % bin size for energies selected according to Rice rule
[Cjo,x2jo]=hist(Muon_Energyjo,bin1jo);
C1jo=Cjo/trapz(x2jo,Cjo);                       % normalization of the histogram
bin8jo=2*(numel(theta_corrjo)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8jo,x8jo]=hist(theta_corrjo,bin8jo);
C9jo=C8jo/trapz(x8jo,C8jo);                     % normalzation of the histogram
bin10jo=2*(numel(Muon_Thetajo)^(1/3));      % bin size for angles selected according to Rice rule
[C10jo,x10jo]=hist(Muon_Thetajo,bin10jo);
C11jo=C10jo/trapz(x10jo,C10jo); 

%%%%%%%%%%%%%%%% tsuji 60
Emints60=3;
Emaxts60=250;
Theta_mints60=59;
Theta_maxts60=61;  

theta1ts60=Theta_mints60*pi/180; theta2ts60=Theta_maxts60*pi/180;             % change from degrees to radians
for i=1:N
tmts60=theta1ts60+(theta2ts60-theta1ts60).*rand();
tm1ts60=2*tmts60/(pi);
theta_corrts60(i)=acos((1-tm1ts60)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2ts60=[theta1ts60:0.01:theta2ts60];                                % zenith angle lower and upper limits
theta3ts60=cos(tm2ts60).^2;                                      % angular distribution modelled as cosine squared
C4ts60=max(theta3ts60);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1ts60=(theta1ts60+(theta2ts60-theta1ts60)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    uts60=rand();                                            % select a random number from uniform U(0,1)
    f1ts60=cos(Theta_m1ts60)^2;                                  % calculate random variable X             
    f2ts60=C4ts60;                                            
    f3ts60=f1ts60/f2ts60;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uts60<=f3ts60                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Thetats60(i)=Theta_m1ts60;                          % accept zenith angle
        break
    end
end

Emts60=[Emints60:0.1:Emaxts60];                                                 % energy lower and upper limits
Epts60=(Emts60+alpha*y0*(sec(Muon_Thetats60(i))-0.1))*(1/rho);
Pmts60=(0.1*cos(Muon_Thetats60(i)).*(1-(alpha*(y0*sec(Muon_Thetats60(i))-100)./(rho.*Epts60)))).^(Bm./((rho.*Epts60+100*alpha)*cos(Muon_Thetats60(i))));
fts60=Alpha.*Pmts60.*(Epts60.^(-kappa))*lambda*bp*jp./(Epts60.*cos(Muon_Thetats60(i))+bp*jp);
C6ts60=max(fts60);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1ts60=(Emints60+(Emaxts60-Emints60)*rand());                           % select a random energy from uniform U(Emin,Emax)
    uts60=rand();                                                % select a random number from uniform U(0,1)
    Ep1ts60=(Em1ts60+alpha*y0*(sec(Muon_Thetats60(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1ts60=(0.1*cos(Muon_Thetats60(i))*(1-(alpha*(y0*sec(Muon_Thetats60(i))-100)/(rho*Ep1ts60))))^(Bm/((rho*Ep1ts60+100*alpha)*cos(Muon_Thetats60(i))));   % calculate probability
    f4ts60=Alpha*Pm1ts60*(Ep1ts60^(-kappa))*lambda*bp*jp/(Ep1ts60*cos(Muon_Thetats60(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5ts60=C6ts60;                                            
    f6ts60=f4ts60/f5ts60;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uts60<=f6ts60                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energyts60(i)=Em1ts60;                                   % accept energy
        break
    end
end
end

Muon_tablets60=[Muon_Thetats60;theta_corrts60;Muon_Energyts60]';

bin1ts60=2*(numel(Muon_Energyts60)^(1/3));      % bin size for energies selected according to Rice rule
[Cts60,x2ts60]=hist(Muon_Energyts60,bin1ts60);
C1ts60=Cts60/trapz(x2ts60,Cts60);                       % normalization of the histogram
bin8ts60=2*(numel(theta_corrts60)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8ts60,x8ts60]=hist(theta_corrts60,bin8ts60);
C9ts60=C8ts60/trapz(x8ts60,C8ts60);                     % normalzation of the histogram
bin10ts60=2*(numel(Muon_Thetats60)^(1/3));      % bin size for angles selected according to Rice rule
[C10ts60,x10ts60]=hist(Muon_Thetats60,bin10ts60);
C11ts60=C10ts60/trapz(x10ts60,C10ts60);  

%%% Datos

kelloggx=[55.6, 70.6, 89.0, 111.3, 140.6, 178.1, 222.6, 309.9, 1067];
kelloggy75=[9.98e-3, 5.41e-3, 3.07e-3, 1.74e-3, 9.50e-4, 4.60e-4, 2.52e-4, 1.09e-4, 2.99e-6];
kelloggy30=[1.58e-2, 7.49e-3, 3.64e-3, 1.80e-3, 9.31e-4, 4.50e-4, 2.25e-4, 7.09e-5, 2.04e-6];


tsujix30=[2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujiy30=[9.57e-4, 6.85e-4, 4.95e-4, 3.70e-4, 2.74e-4, 2.06e-4, 1.61e-4, 1.37e-4, 7.74e-5, 3.61e-5, 1.88e-5, 1.17e-5, 5.72e-6, 2.79e-6, 1.57e-6, 8.42e-7, 5.97e-7, 5.20e-7, 2.16e-7, 8.10e-8, 4.29e-8, 1.31e-8];
tsujix60=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujiy60=[2.22e-4, 1.91e-4, 1.42e-4, 1.07e-4, 9.94e-5, 8.01e-5, 6.35e-5, 4.71e-5, 2.35e-5, 1.18e-5, 7.46e-6, 4.08e-6, 1.78e-6, 1.64e-6, 1.08e-6, 2.19e-7, 4.26e-7, 2.38e-7, 7.37e-8, 1.17e-7, 3.49e-8];

jokix75=[1.12, 1.41, 1.76, 2.29, 2.85, 3.57, 4.49, 5.66, 7.12, 8.95, 11.26, 14.17, 17.83, 22.44, 28.23, 35.52, 44.70, 56.25, 70.79, 89.09, 112.13, 141.12, 177.62, 223.56, 281.39, 354.18, 445.81, 625.46, 990.22];
jokiy75=[2.97e-5,  3.11e-5, 3.263e-5, 3.168e-5, 3.068e-5, 2.848e-5, 2.606e-5, 2.332e-5, 2.021e-5, 1.705e-5, 1.373e-5, 1.070e-5, 8.029e-6, 5.708e-6, 3.964e-6, 2.637e-6, 1.712e-6, 1.064e-6, 6.479e-7, 3.751e-7, 2.208e-7, 1.234e-7, 6.632e-8, 3.735e-8, 1.869e-8, 9.689e-9, 5.436e-9, 1.770e-9, 3.982e-10];




%%%% graficas
close all;
figure(1)
hold on
loglog(x2,C1/5,'^k')
loglog(x2kell,C1kell/5,'^r')
loglog(x2ts30,C1ts30/1000,'^g')
loglog(x2jo,C1jo/1000,'^b')
loglog(x2ts60,C1ts60/1000,'^m')
hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

hold on
plot(kelloggx, kelloggy30,'*k', linewidth=1)
plot(kelloggx, kelloggy75,'*r', linewidth=1)
plot(tsujix30, tsujiy30,'*g', linewidth=1)
plot(jokix75, jokiy75,'*b', linewidth=1)
plot(tsujix60, tsujiy60,'*m', linewidth=1)
hold off
legend({'Gen Kellogg (30°)','Gen Kellogg (75°)','Gen Tsuji (30°)','Gen Jokisch (75°)','Gen Tsuji (60°)', 'Kellogg (30°)', 'Kellogg (75°)','Tsuji (30°)','Jokisch (75°)','Tsuji (60°)'},Location="best")
xlabel('Momento del muón [GeV/c]')
ylabel('Intensidad diferencial [cm^{-2}sr^{-1}sec^{-1}(GeV/c)^{-1}]')
