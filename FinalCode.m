% functionSynapseandflux
%Parameters
L=0.02; %µm
tfinal=1e-7; %seconds
D=400; %µm^2/s
tpoints=501;
xpoints=41;
%m=0 denotes that our problem is in a slab
m = 0;
%Create time and distance vectors
x = linspace(0,L,xpoints);
t = linspace(0,tfinal,tpoints);
%Call pdepe to solve the equation
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
%Solution is first component of sol
C = sol(:,:,1);
%Estimating flux at boundary
Flux=zeros(1,length(t));
for i=2:501
dFlux=-D*((C(i,41)-C(i,40))/(L/xpoints))*(tfinal/tpoints);
Flux(i)=Flux(i-1)+dFlux;
end
%Plotting Flux vs time
figure
plot(t,Flux*6.022e23)
xlabel('Time (seconds)')
ylabel('Bound Receptor Surface Concentration (molecules/µm^2)')
title('ACh Bound at Post-Synaptic Membrane')
%Create 3d plot
figure
surf(x,t,C,'Edgecolor','None')
title('AChConc in the Synaptic Cleft')
xlabel('Distance (µm)')
ylabel('Time (seconds)')
zlabel('Concentration (mol/m^3)')

%PDE function (define our PDE)
function [c,f,s] = pdex1pde(x,t,Conc,DCDx)
15
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
Cs=6e-13; %mol/µm^2, concentration pulse into synapse
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
beta=4.7e22; %µm^3/mol*s Receptor binding rate
%dCdX=0 at x=0
pl = 0;
ql = 1;
%-beta*Conc-D*dCdx=0 at x=L (Robin BC)
pr = -beta*ur;
qr = -1;
end