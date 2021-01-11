% NOTE: Assume h-bar=m=1
% Exact solution for particle in 1D box
clc, clear
nx=input('Quantum number n = ');    % Quantum number along x-axis
x=0:0.01:1;                         % Spatial domain
t=0:0.001:0.5;                      % Temporal domain
psi=zeros(1,length(x));

title(sprintf('Exact solution for a particle in an infinite potential well with quantum number n=%d',nx));
axis([0 1 -2 2]);
xlabel('Length of 1D box');
hold on; grid on;
Enx=0.5*(nx*pi)^2;
for i=1:length(t)
    
    % Plot for particle time-dependent wave function
    [Re,Im]=Psi(nx,Enx,t(i),x);
    box=plot(x,Re,'r');
    
    % Plot of energy density
    rho=Re.^2-Im.^2;
    e=plot(x,rho,'b');
    
    legend('Wave function (\Psi)','Energy density (\rho)');
    txt=sprintf('At t = %.3f',t(i));
    time=text(0.02,1.85,txt);
    drawnow
    if i~=length(t)
        delete(box); delete(time); delete(e);
    end
end

%% List of functions requried for above
function [Re,Im]=Psi(nx,Enx,time,x) % Wave-function exact solver
[cosine, sine]=euler2polar(-Enx*time);
Re=sqrt(2)*cosine.*sin(nx*pi*x);
Im=sqrt(2)*sine.*sin(nx*pi*x);
end

% Convert Euler equation to polar
function [cosine, sine]=euler2polar(theta)
cosine=cos(theta);
if theta<0
    sine=-1i*sin(abs(theta));
else
    sine=1i*sin(abs(theta));
end
end

