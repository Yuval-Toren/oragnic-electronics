% Numerical solution (FDM)
clc, clear
dx=0.01;                  % Spatial discretization
x=0:dx:1;                 % Spatial domain
C=0.01;                   % User-defined value until stability is reached
dt=C*dx^2;                % Time discretization
t=0:dt:0.5;               % Temporal domain
dpsi=zeros(1,length(x));  % Derivative of psi w.r.t. time
% Initial conditions - where t(1)=0
nx=input('Quantum number n = '); 
Enx=0.5*(nx*pi)^2;
[Re2,Im2]=Psi(nx,Enx,t(1),x);
psi=Re2+Im2;

% Main solver based on FDM & RK3 method
axis([0 1 -2 2]);
xlabel('Length of 1D box');
hold on; grid on;
axis manual;
for i=1:length(t)
    % B.C.s at ends are zero cause infinite potential
    psi(1)=0; psi(length(x))=0; % Dirichlet boundary conditions
                        
    dpsi=dpsidt(psi,dx);
    k1=dt*dpsi;
    
    dpsi2=dpsidt(psi+k1/2,dx);
    k2=dt*dpsi2;
    
    dpsi3=dpsidt(psi+k2/2,dx);
    k3=dt*dpsi3;
        
    psi=psi+(1/6)*k1+(2/3)*k2+(1/6)*k3;
    
    if max(abs(psi)) > 1.5   % Check for divergence
        sprintf('Error!!! At i=%d',i)
        return
    end
    
    % Simulation
    if rem(i,10/dx)==0 || i==1
        box2=plot(x,real(psi),'r');
        
        % Plot of energy density
        rho=real(psi).^2+imag(psi).^2;
        e=plot(x,rho,'b');
    
        legend('Wave function (\Psi)','Energy density (\rho)');
        txt=sprintf('At t = %.3f',t(i));
        time=text(0.02,1.85,txt);
        drawnow;
        if length(t)/1000-(i/1000)>1/1000
            delete(box2); delete(time); delete(e);
        end
    end
end

function [dpsi]=dpsidt(psi,dx)
len=length(psi);
dpsi=zeros(1,len);
for j=2:len-1
    dpsi(j)=0.5i*((psi(j+1)-2*psi(j)+psi(j-1))/dx^2);
end
% Boundary conditions (differential)
dpsi(1)=0.5i*((psi(3)-2*psi(2)+psi(1))/dx^2);             % For left boundary
dpsi(len)=0.5i*((psi(len)-2*psi(len-1)+psi(len-2))/dx^2); % For right boundary
end

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
