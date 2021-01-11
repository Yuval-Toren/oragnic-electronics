% Exact solution
clc, clear
x=-10:0.1:10;                       % Spatial domain
t=0:0.01:100;                     % Temporal domain
psi=zeros(1,length(x));
H=1;                              % Height of vertical axis for plotting

% Constants for Hermite polynomial 
m=1; hbar=1;                      % NOTE: Assume h-bar=m=1
n=input('Energy level, n = ');    % Energy level
omega=1;
k=m*omega^2;
E=hbar*omega*(n+0.5);

title(sprintf('Exact solution at energy level n=%d',n));
axis([x(1) x(length(x)) -H H]);
xlabel('Length');
hold on; grid on;

% Main solver
Enx=0.5*(n*pi)^2;
sprintf('Note that the Hermite polynomial used for this is %d-order',n)
for i=1:length(t)
    
    % Plot for particle time-dependent wave function
    [Re,Im]=Psi(k,E,t(i),x,n);
    box=plot(x,Re,'r');
    
    % Plot of probability density
    rho=(Re+Im).*conj(Re+Im);
    e=plot(x,rho,'b');
    
    legend('Wave function (\Psi)','Probability density (|\Psi|^2)');
    txt=sprintf('At t = %.2f',t(i));
    time=text(0.05+x(1),H-0.15,txt);
    drawnow
    if i~=length(t)
        delete(box); delete(time); delete(e);
    end
end

%% List of functions required for above
function [Re,Im]=Psi(k,E,time,x,n) % Wave-function exact solver
% Hermite polynomial
syms k
H=factorial(n).*symsum( (((-1)^(k))./(factorial(k)*factorial(n-2*k))).*(2*x).^(n-2*k),k,0,floor(n/2) );

[cosine, sine]=euler2polar(-E*time);
Re=sqrt(1/(sqrt(pi)*2^n*factorial(n)))*cosine.*exp(-0.5*x.^2).*H;
Im=sqrt(1/(sqrt(pi)*2^n*factorial(n)))*sine.*exp(-0.5*x.^2).*H;
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

