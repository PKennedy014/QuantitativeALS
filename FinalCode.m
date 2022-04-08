clc,clear,
% BME4409 Final Project
% Dr. Ferrall-Fairbanks
% Members: Jeffrey Chen, Vanessea Cruz, Estefania Hernandez, Patrick Kennedy
%Modeling the effect of ALS on muscle contraction

%Parameters
L=40*10^-3; %Cleft distance µm 
tfinal=2e-7; %seconds. For graph
D=400; %ACh difusivity in NMJ µm^2/s
tpoints=501;
xpoints=41;

%m=0 denotes that our problem is in a slab
m = 0;

%Create time and distance vectors
x = linspace(0,L,xpoints);
t = linspace(0,tfinal,tpoints);

%Maximum bound receptor surface concentration
BRSCmax = zeros(10,1);
for j = 1:10
    global ALS
    ALS = j; %y being severity of ALS scale (1-10); with 1 being normal

    %Call pdepe to solve the equation
    sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
    %Solution is first component of sol
    C = sol(:,:,1);
    %Estimating flux at boundary
    Flux=zeros(1,length(t));
    for i=2:tpoints
    dFlux=-D*((C(i,41)-C(i,40))/(L/xpoints))*(tfinal/tpoints);
    Flux(i)=Flux(i-1)+dFlux;
    end
    %Bound Receptor surface concentration
    BRSC = Flux*6.022e23;
    BRSCmax(j) = max(BRSC);
    
%     %Plotting BRSC vs time
%     figure
%     plot(t,BRSC)
%     xlabel('Time (seconds)')
%     ylabel('Bound Receptor Surface Concentration (molecules/µm^2)')
%     title('ACh Bound at Post-Synaptic Membrane')
end

%Plotting BRSCmax vs ALS
subplot(1,2,1)
bar(1:10, BRSCmax)
xlabel('Severity of ALS (scale 1-10)')
ylabel('Max Bound Receptor Surface Concentration (molecules/µm^2)')
title('ALS Effect on Bound Nicotinic Receptors on Motor End Plate')

%Plotting BRSCmax vs. % muscle contraction
BRSCmaxSmooth = (0:0.01:max(BRSCmax));
MusCon = 100./(1+exp(-(BRSCmaxSmooth-35)./4));
subplot(1,2,2)
plot(BRSCmaxSmooth,MusCon)
xlabel('Max Bound Receptor Surface Concentration (molecules/µm^2)')
ylabel('% Muscle contraction')
title('Bound Nicotinic Receptors on Motor End Plate vs. %Contraction')

%Plotting % muscle contraction vs. ALS
figure; MusCon2 = 100./(1+exp(-(BRSCmax-35)./4));
plot(1:10,MusCon2)
title('ALS Effect on Muscle Contraction')
xlabel('Severity of ALS (scale 1-10)')
ylabel('% Muscle contraction')

%PDE function (define our PDE)
function [c,f,s] = pdex1pde(x,t,Conc,DCDx)
D=400; %µm^2/s
k=23e7; %1/s, kcat of AChE
%Time derivative coefficient
c = 1;
%x derivative coefficient (includes first derivative to make a second
%derivative in x)
f = D*DCDx;
%Forcing function coefficient
G=-k*Conc; %mol/µm^3 s
s = G;
end

%IC function
function C0 = pdex1ic(x)
L=0.02;
global ALS
Cs=(1/ALS)*6e-6; %mol/µm^2, concentration pulse into synapse
%Delta function initial condition, concentration is Cs at x=0 and 0
%elsewhere
if x==0
 C0=Cs;
else
 C0=0;
end
end

%BC function
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
D=400; %µm^2/s;
global ALS
beta=(1/ALS)*4.7e22; %µm^3/mol*s Receptor binding rate
%dCdX=0 at x=0
pl = 0;
ql = 1;
%-beta*Conc-D*dCdx=0 at x=L (Robin BC)
pr = -beta*ur;
qr = -1;
end