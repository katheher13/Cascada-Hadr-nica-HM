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

%%% Nandi:
Emin=5;
Emax=1200;
Theta_min=0;
Theta_max=0.3;  
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

%%%%%%%%%%%%%%%%%%%%% DATOS
baberx = [3.55, 3.75, 4.17, 4.71, 5.26, 6.01, 7.57, 9.30, 11.60, 15.22, 19.20, 24.00, 33.5, 50.0, 81.0, 127.0, 266.0, 810.0];
babery = [6.27e-4, 5.84e-4, 5.03e-4, 4.20e-4, 3.53e-4, 2.82e-4, 1.86e-4, 1.24e-4, 7.77e-5, 4.22e-5, 2.42e-5, 1.39e-5, 5.78e-6, 1.90e-6, 4.59e-7, 1.14e-7, 1.00e-8, 2.11e-10];

nandix=[5,10,15,20,40, 60, 80, 100, 150, 200, 400, 600, 800, 1000 ];
nandiy=[4.70e-4, 1.30e-4, 5.39e-5, 2.73e-5, 4.55e-6, 1.46e-6, 6.24e-7, 3.18e-7, 8.97e-8, 3.56e-8, 3.53e-9, 8.75e-10, 3.25e-10, 1.48e-10];

allk1x=[0.34, 0.64, 1.32, 2.86, 4.91, 1.11, 0.29, 0.50, 0.79, 1.24, 1.95, 3.17, 4.91, 7.76, 11.4, 14.8, 20.5, 31.4, 52.3, 93.0, 175.0, 329.0, 642.0];
allk1y=[3.92e-3, 3.59e-3, 2.57e-3, 1.11e-3, 4.13e-4, 2.90e-3, 3.57e-3, 3.70e-3, 3.41e-3, 2.73e-3, 1.73e-3, 7.92e-4, 4.24e-4, 1.84e-4, 1.13e-4, 6.04e-5, 2.51e-5, 8.01e-6, 1.89e-6, 3.38e-7, 5.19e-8, 7.84e-9, 6.40e-10];


%length(allk1y)
%length(allk1x)
%%%%%%%%%%% GRAFICAS
close all;
figure(1)
hold on
loglog(x2,C1/1000,'^b')
loglog(x2baber,C1baber/1000,'^g')
loglog(x2All,C1All/1000,'^m')
hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

hold on
plot(nandix, nandiy,'*b', linewidth=1)
plot(baberx, babery,'*g', linewidth=1)
plot(allk1x, allk1y,'*m', linewidth=1)
hold off
legend({'Gen Nandi','Gen Baber','Gen Allkofer','Nandi','Baber', 'Allkofer'},Location="best")
xlabel('Momento del muón [GeV/c]')
ylabel('Intensidad diferencial [cm^{-2}sr^{-1}sec^{-1}(GeV/c)^{-1}]')