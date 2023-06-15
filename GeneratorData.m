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

%%%%  BABER 
Eminbaber=3;
Emaxbaber=1000;
Theta_minbaber=0;
Theta_maxbaber=0;  

theta1baber=Theta_minbaber*pi/180; theta2baber=Theta_maxbaber*pi/180;             % change from degrees to radians
for i=1:N
tmbaber=theta1baber+(theta2baber-theta1baber).*rand();
tm1baber=2*tmbaber/(pi);
theta_corrbaber(i)=acos((1-tm1baber)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2baber=[theta1baber:0.01:theta2baber];                                % zenith angle lower and upper limits
theta3baber=cos(tm2baber).^2;                                      % angular distribution modelled as cosine squared
C4baber=max(theta3baber);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1baber=(theta1baber+(theta2baber-theta1baber)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    ubaber=rand();                                            % select a random number from uniform U(0,1)
    f1baber=cos(Theta_m1baber)^2;                                  % calculate random variable X             
    f2baber=C4baber;                                            
    f3baber=f1baber/f2baber;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ubaber<=f3baber                                             % if u<=f(y)/g(y) then set X=Y
        Muon_Thetababer(i)=Theta_m1baber;                          % accept zenith angle
        break
    end
end

Embaber=[Eminbaber:0.1:Emaxbaber];                                                 % energy lower and upper limits
Epbaber=(Embaber+alpha*y0*(sec(Muon_Thetababer(i))-0.1))*(1/rho);
Pmbaber=(0.1*cos(Muon_Thetababer(i)).*(1-(alpha*(y0*sec(Muon_Thetababer(i))-100)./(rho.*Epbaber)))).^(Bm./((rho.*Epbaber+100*alpha)*cos(Muon_Thetababer(i))));
fbaber=Alpha.*Pmbaber.*(Epbaber.^(-kappa))*lambda*bp*jp./(Epbaber.*cos(Muon_Thetababer(i))+bp*jp);
C6baber=max(fbaber);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1baber=(Eminbaber+(Emaxbaber-Eminbaber)*rand());                           % select a random energy from uniform U(Emin,Emax)
    ubaber=rand();                                                % select a random number from uniform U(0,1)
    Ep1baber=(Em1baber+alpha*y0*(sec(Muon_Thetababer(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1baber=(0.1*cos(Muon_Thetababer(i))*(1-(alpha*(y0*sec(Muon_Thetababer(i))-100)/(rho*Ep1baber))))^(Bm/((rho*Ep1baber+100*alpha)*cos(Muon_Thetababer(i))));   % calculate probability
    f4baber=Alpha*Pm1baber*(Ep1baber^(-kappa))*lambda*bp*jp/(Ep1baber*cos(Muon_Thetababer(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5baber=C6baber;                                            
    f6baber=f4baber/f5baber;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if ubaber<=f6baber                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_Energybaber(i)=Em1baber;                                   % accept energy
        break
    end
end
end

Muon_tablebaber=[Muon_Thetababer;theta_corrbaber;Muon_Energybaber]';

bin1baber=2*(numel(Muon_Energybaber)^(1/3));      % bin size for energies selected according to Rice rule
[Cbaber,x2baber]=hist(Muon_Energybaber,bin1baber);
C1baber=Cbaber/trapz(x2baber,Cbaber);                       % normalization of the histogram
bin8baber=2*(numel(theta_corrbaber)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8baber,x8baber]=hist(theta_corrbaber,bin8baber);
C9baber=C8baber/trapz(x8baber,C8baber);                     % normalzation of the histogram
bin10baber=2*(numel(Muon_Thetababer)^(1/3));      % bin size for angles selected according to Rice rule
[C10baber,x10baber]=hist(Muon_Thetababer,bin10baber);
C11baber=C10baber/trapz(x10baber,C10baber);  
%%%%%%%%%%%%%%%%%%%%% ALLKOFER
EminAll=0.2;
EmaxAll=1000;
Theta_minAll=0;
Theta_maxAll=0;  

theta1All=Theta_minAll*pi/180; theta2All=Theta_maxAll*pi/180;             % change from degrees to radians
for i=1:N
tmAll=theta1All+(theta2All-theta1All).*rand();
tm1All=2*tmAll/(pi);
theta_corrAll(i)=acos((1-tm1All)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform

tm2All=[theta1All:0.01:theta2All];                                % zenith angle lower and upper limits
theta3All=cos(tm2All).^2;                                      % angular distribution modelled as cosine squared
C4All=max(theta3All);                                          % maximum value of angular distribution

for k=1:1000000                       
    Theta_m1All=(theta1All+(theta2All-theta1All)*rand());            % select a random energy from uniform U(Theta_min,Theta_max)
    uAll=rand();                                            % select a random number from uniform U(0,1)
    f1All=cos(Theta_m1All)^2;                                  % calculate random variable X             
    f2All=C4All;                                            
    f3All=f1All/f2All;                                            % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uAll<=f3All                                             % if u<=f(y)/g(y) then set X=Y
        Muon_ThetaAll(i)=Theta_m1All;                          % accept zenith angle
        break
    end
end

EmAll=[EminAll:0.1:EmaxAll];                                                 % energy lower and upper limits
EpAll=(EmAll+alpha*y0*(sec(Muon_ThetaAll(i))-0.1))*(1/rho);
PmAll=(0.1*cos(Muon_ThetaAll(i)).*(1-(alpha*(y0*sec(Muon_ThetaAll(i))-100)./(rho.*EpAll)))).^(Bm./((rho.*EpAll+100*alpha)*cos(Muon_ThetaAll(i))));
fAll=Alpha.*PmAll.*(EpAll.^(-kappa))*lambda*bp*jp./(EpAll.*cos(Muon_ThetaAll(i))+bp*jp);
C6All=max(fAll);                                                           % maximum value of phenomenological model

for m=1:1000000                       
    Em1All=(EminAll+(EmaxAll-EminAll)*rand());                           % select a random energy from uniform U(Emin,Emax)
    uAll=rand();                                                % select a random number from uniform U(0,1)
    Ep1All=(Em1All+alpha*y0*(sec(Muon_ThetaAll(i))-0.1))*(1/rho);     % calculate pion energy
    Pm1All=(0.1*cos(Muon_ThetaAll(i))*(1-(alpha*(y0*sec(Muon_ThetaAll(i))-100)/(rho*Ep1All))))^(Bm/((rho*Ep1All+100*alpha)*cos(Muon_ThetaAll(i))));   % calculate probability
    f4All=Alpha*Pm1All*(Ep1All^(-kappa))*lambda*bp*jp/(Ep1All*cos(Muon_ThetaAll(i))+bp*jp);    % calculate intensity based on sampled energy and angle             
    f5All=C6All;                                            
    f6All=f4All/f5All;                                                 % f(y)/g(y)=f1/f2 where g(y)>f(y)
    if uAll<=f6All                                                  % if u<=f(y)/g(y) then set X=Y
        Muon_EnergyAll(i)=Em1All;                                   % accept energy
        break
    end
end
end

Muon_tableAll=[Muon_ThetaAll;theta_corrAll;Muon_EnergyAll]';

bin1All=2*(numel(Muon_EnergyAll)^(1/3));      % bin size for energies selected according to Rice rule
[CAll,x2All]=hist(Muon_EnergyAll,bin1All);
C1All=CAll/trapz(x2All,CAll);                       % normalization of the histogram
bin8All=2*(numel(theta_corrAll)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8All,x8All]=hist(theta_corrAll,bin8All);
C9All=C8All/trapz(x8All,C8All);                     % normalzation of the histogram
bin10All=2*(numel(Muon_ThetaAll)^(1/3));      % bin size for angles selected according to Rice rule
[C10All,x10All]=hist(Muon_ThetaAll,bin10All);
C11All=C10All/trapz(x10All,C10All);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kellogg

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
Theta_minjo=62;
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

%%%%%%%%%%%%%%%%%%%%% DATOS
baberx = [3.55, 3.75, 4.17, 4.71, 5.26, 6.01, 7.57, 9.30, 11.60, 15.22, 19.20, 24.00, 33.5, 50.0, 81.0, 127.0, 266.0, 810.0];
babery = [6.27e-4, 5.84e-4, 5.03e-4, 4.20e-4, 3.53e-4, 2.82e-4, 1.86e-4, 1.24e-4, 7.77e-5, 4.22e-5, 2.42e-5, 1.39e-5, 5.78e-6, 1.90e-6, 4.59e-7, 1.14e-7, 1.00e-8, 2.11e-10];

nandix=[5,10,15,20,40, 60, 80, 100, 150, 200, 400, 600, 800, 1000 ];
nandiy=[4.70e-4, 1.30e-4, 5.39e-5, 2.73e-5, 4.55e-6, 1.46e-6, 6.24e-7, 3.18e-7, 8.97e-8, 3.56e-8, 3.53e-9, 8.75e-10, 3.25e-10, 1.48e-10];

allk1x=[0.34, 0.64, 1.32, 2.86, 4.91, 1.11, 0.29, 0.50, 0.79, 1.24, 1.95, 3.17, 4.91, 7.76, 11.4, 14.8, 20.5, 31.4, 52.3, 93.0, 175.0, 329.0, 642.0];
allk1y=[3.92e-3, 3.59e-3, 2.57e-3, 1.11e-3, 4.13e-4, 2.90e-3, 3.57e-3, 3.70e-3, 3.41e-3, 2.73e-3, 1.73e-3, 7.92e-4, 4.24e-4, 1.84e-4, 1.13e-4, 6.04e-5, 2.51e-5, 8.01e-6, 1.89e-6, 3.38e-7, 5.19e-8, 7.84e-9, 6.40e-10];

ayrex=[20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500];
ayrey=[2.43e-5, 8.80e-6, 4.06e-6, 2.20e-6, 1.31e-6, 8.43e-7, 5.71e-7, 4.04e-7, 2.96e-7, 8.68e-8, 3.56e-8, 9.77e-9, 3.80e-9, 1.79e-9];

kelloggx=[55.6, 70.6, 89.0, 111.3, 140.6, 178.1, 222.6, 309.9, 1067];
kelloggy75=[9.98e-3, 5.41e-3, 3.07e-3, 1.74e-3, 9.50e-4, 4.60e-4, 2.52e-4, 1.09e-4, 2.99e-6];
kelloggy30=[1.58e-2, 7.49e-3, 3.64e-3, 1.80e-3, 9.31e-4, 4.50e-4, 2.25e-4, 7.09e-5, 2.04e-6];

%tsujix30=[1.8, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujix75=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsuji75=[4.59e-5, 4.16e-5, 3.67e-5, 3.33e-5, 2.97e-5, 2.76e-5, 2.43e-5, 1.98e-5, 1.09e-5, 7.14e-6, 4.82e-6, 2.51e-6, 1.50e-6, 7.19e-7, 4.77e-7, 2.79e-7, 2.43e-7, 2.26e-7, 8.98e-8, 1.27e-8, 4.43e-8];
%tsuji30=[9.57e-4, 6.85e-4, 4.95e-4, 3.70e-4, 2.74e-4, 2.06e-4, 1.61e-4, 1.37e-4, 7.74e-5, 3.61e-5, 1.88e-5, 1.17e-5, 5.72e-6, 2.79e-6, 1.57e-6, 8.42e-7, 5.97e-7, 5.20e-7, 2.16e-7, 8.10e-8, 4.29e-8, 1.31e-8];
tsujix0=[1.8, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujiy0=[1.51e-3, 1.21e-3, 7.93e-4, 5.73e-4, 3.93e-4, 3.01e-4, 2.30e-4, 1.78e-4, 1.38e-4, 8.37e-5, 3.75e-5, 2.04e-5, 1.17e-5, 6.12e-6, 3.21e-6, 1.39e-6, 8.68e-7, 8.07e-7, 3.73e-7, 1.69e-7, 1.00e-7, 2.84e-8, 3.59e-8];
tsujix30=[2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujiy30=[9.57e-4, 6.85e-4, 4.95e-4, 3.70e-4, 2.74e-4, 2.06e-4, 1.61e-4, 1.37e-4, 7.74e-5, 3.61e-5, 1.88e-5, 1.17e-5, 5.72e-6, 2.79e-6, 1.57e-6, 8.42e-7, 5.97e-7, 5.20e-7, 2.16e-7, 8.10e-8, 4.29e-8, 1.31e-8];
tsujix60=[3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 35, 45, 55, 65, 75, 85, 95, 125, 175, 225];
tsujiy60=[2.22e-4, 1.91e-4, 1.42e-4, 1.07e-4, 9.94e-5, 8.01e-5, 6.35e-5, 4.71e-5, 2.35e-5, 1.18e-5, 7.46e-6, 4.08e-6, 1.78e-6, 1.64e-6, 1.08e-6, 2.19e-7, 4.26e-7, 2.38e-7, 7.37e-8, 1.17e-7, 3.49e-8];


jokix75=[1.12, 1.41, 1.76, 2.29, 2.85, 3.57, 4.49, 5.66, 7.12, 8.95, 11.26, 14.17, 17.83, 22.44, 28.23, 35.52, 44.70, 56.25, 70.79, 89.09, 112.13, 141.12, 177.62, 223.56, 281.39, 354.18, 445.81, 625.46, 990.22];
jokiy75=[2.97e-5,  3.11e-5, 3.263e-5, 3.168e-5, 3.068e-5, 2.848e-5, 2.606e-5, 2.332e-5, 2.021e-5, 1.705e-5, 1.373e-5, 1.070e-5, 8.029e-6, 5.708e-6, 3.964e-6, 2.637e-6, 1.712e-6, 1.064e-6, 6.479e-7, 3.751e-7, 2.208e-7, 1.234e-7, 6.632e-8, 3.735e-8, 1.869e-8, 9.689e-9, 5.436e-9, 1.770e-9, 3.982e-10];

%length(allk1y)
%length(allk1x)
%%%%%%%%%%% GRAFICAS
close all;
figure(1)
hold on
loglog(x2baber,C1baber/1000,'^k')
loglog(x2All,C1All/1000,'^r')
%loglog(x2kell,C1kell/5,':g')
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
plot(baberx, babery,'*k', linewidth=1)
plot(allk1x, allk1y,'*r', linewidth=1)
%plot(kelloggx, kelloggy75,'--g', linewidth=0.5)
plot(tsujix30, tsujiy30,'*g', linewidth=1)
plot(jokix75, jokiy75,'*b', linewidth=1)
plot(tsujix60, tsujiy60,'*m', linewidth=1)
hold off
legend({'Gen Baber (0°)','Gen Allkofer (0°)','Gen Tsuji (30°)','Gen Jokisch (75°)','Gen Tsuji (60°)', 'Baber (0°)', 'Allkofer (0°)','Tsuji (30°)','Jokisch (75°)','Tsuji (60°)'},Location="best")
xlabel('Momento del muón [GeV/c]')
ylabel('Intensidad diferencial [cm^{-2}sr^{-1}sec^{-1}(GeV/c)^{-1}]')