%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Muon generator matlab file (Version 2).
%   
%   This file provides a cosmic ray muon sampling capability for use in 
%   Monte Carlo simulations. The muon energy distribution is described by 
%   the Smith & Duller (1959) phenomenological model that captures the main 
%   characteristics of the experimentally measured spectrum. 
%   Statistical algorithms are then employed for generating random samples. 
%   The inverse transform provides a means to generate samples from the 
%   muon angular distribution (corrected for solid angle effect). The Acceptance-Rejection 
%   method is used to generate samples from the actual muon angular
%   distribution (i.e., cosine squared). The Acceptance-Rejection 
%   method is also used to provide the energy component. Final output is a 
%   look-up table with three columns that contains the sampled angles and energies.
%
%   Instructions: Run the file. On the command line set the minimum and 
%   maximum muon energy, the minimum and maximum zenith angle and the 
%   number of muons to be sampled. The output "Muon_table" matrix contains 
%   the angles and energies of the sampled muons.
%
%   This file is free for use. More details, examples and validation results 
%   can be found in the journal paper:
%   
%   S. Chatzidakis, S. Chrysikopoulou, L.H. Tsoukalas (2015)
%   "Developing a cosmic ray muon sampling capability for muon tomography
%   and monitoring applications", Nuclear Instruments and Methods in 
%   Phyics Research Section A, Vol. 804, pp. 33-42.
%
%   S. Chatzidakis, L.H. Tsoukalas (2015)
%   "A Geant4-MATLAB muon generator for Monte-Carlo simulations", 
%   URL: https://engineering.purdue.edu/~aisl/Stylianos_Publications.html
%
%   Users are kindly requested to cite the above journal paper in their
%   publications. Thank you!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Set minimum and maximum muon energy (GeV) (Range 1-100 GeV)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Emin=input('Enter minimum muon energy (GeV):');
Emax=input('Enter maximum muon energy (GeV):');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Set minimum and maximum muon zenith angle (Range 0-90o)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta_min=input('Enter minimum muon zenith angle (degrees):');
Theta_max=input('Enter maximum muon zenith angle (degrees):');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Set the desirable number of muons to be sampled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=input('Enter the number of muons:');                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Phenomenological model constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alpha=0.002382;              % constant A
lambda=120;                  % absorption mean free path 120 g/cm2
kappa=2.645;                 % exponent (-)
bp=0.771; jp=148.16;         % correction factor (-); factor (GeV)
alpha = 0.0025;              % muon energy loss in GeV/g/cm2
rho = 0.76;                  % fraction of pion energy that is transferred to muon
y0=1000;                     % atmoshperic depth g/cm2
bm=0.8;Bm=1.041231831;       % correction factor (-); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inverse trasform and Accept-Reject method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Step 1
%   Inverse trasform to sample from muon zenith 
%   angle - corrected for solid angle effect
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1=Theta_min*pi/180;theta2=Theta_max*pi/180;             % change from degrees to radians
for i=1:N
tm=theta1+(theta2-theta1).*rand();
tm1=2*tm/(pi);
theta_corr(i)=acos((1-tm1)^(1/3));                                % select muon angle (in radians) - corrected for solid angle effect - using inverse transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Step 2
%    Accept-Reject method to sample muon zenith angle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate maximum value of angular distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tm2=[theta1:0.01:theta2];                                % zenith angle lower and upper limits
theta3=cos(tm2).^2;                                      % angular distribution modelled as cosine squared
C4=max(theta3);                                          % maximum value of angular distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Actual angular muon spectrum
%   (i.e., cosine squared)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Step 3
%    Accept-Reject method to sample muon energy 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate maximum value of phenomenological model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Em=[Emin:0.1:Emax];                                                 % energy lower and upper limits
Ep=(Em+alpha*y0*(sec(Muon_Theta(i))-0.1))*(1/rho);
Pm=(0.1*cos(Muon_Theta(i)).*(1-(alpha*(y0*sec(Muon_Theta(i))-100)./(rho.*Ep)))).^(Bm./((rho.*Ep+100*alpha)*cos(Muon_Theta(i))));
f=Alpha.*Pm.*(Ep.^(-kappa))*lambda*bp*jp./(Ep.*cos(Muon_Theta(i))+bp*jp);
C6=max(f);                                                           % maximum value of phenomenological model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Muon energy sampling 
%  (for the previously selected muon angle Muon_Theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Create a lookup table: 1st column is sampled muon angle (radians), 
%   second column is sampled muon angle -corrected for solid angle effect (radians)
%   third colum is sampled muon energy (GeV)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Muon_table=[Muon_Theta;theta_corr;Muon_Energy]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Set bin size and produce histograms for plotting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin1=2*(numel(Muon_Energy)^(1/3));      % bin size for energies selected according to Rice rule
[C,x2]=hist(Muon_Energy,bin1);
C1=C/trapz(x2,C);                       % normalization of the histogram
bin8=2*(numel(theta_corr)^(1/3));            % bin size for angles (corrected for solid angle effect) selected according to Rice rule
[C8,x8]=hist(theta_corr,bin8);
C9=C8/trapz(x8,C8);                     % normalzation of the histogram
bin10=2*(numel(Muon_Theta)^(1/3));      % bin size for angles selected according to Rice rule
[C10,x10]=hist(Muon_Theta,bin10);
C11=C10/trapz(x10,C10);                 % normalzation of the histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot histograms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)                       
plot(x8,C9,'o')                         % plot distribution of zenith angles - corrected for the solid angle effect.
xlabel('Zenith angle - Corrected for solid angle effect (rad)')
ylabel('Arbitrary units')
figure(2)                       
plot(x10,C11,'o')                         % plot distribution of zenith angles.
xlabel('Zenith angle (rad)')
ylabel('Arbitrary units')
figure(3)
loglog(x2,C1,'o')                         % plot distribution of muon energies
xlabel('Energy (GeV)')
ylabel('Arbitrary units')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Clear workspace
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Muon_table    % clear workspace
