%% 2 a
U = -1.922*(10^(-18));%eV
E = -7.718*(10^(-19));%eV
m = 9.10*10^(-31);%Kg
L = 1;%m
hbar = 1.054*10^(-34);%eV*s
lambda = 10^(-10)*sqrt(-2*m*E/(hbar^2));
k = 10^(-10)*sqrt(2*m*(E-U)/(hbar^2));
ee = exp(-lambda*L/2);
ss = sin(k*L/2);
cc = cos(k*L/2);
A1 = 1;
A2 = A1*(1+(tan(k*L/2))^2)/(2*(lambda/k)*tan(k*L/2)+...
    1-(tan(k*L/2))^2);
C = (ee*(A1+A2))/(2*cc);
B = (ee*(A2-A1))/(2*ss);
xi = 8;

fun = @(x) (abs((B*sin(k*(x-xi))+C*cos(k*(x-xi))))).^2;
Sii = integral(fun,7.5,8.5)
%% 1 b 2
% U = -2.035*(10^(-18));%eV
% E = -9.1324e-19;%eV
% m = 9.10*10^(-31);%Kg
% L = 1.1*10^(-10);%m
% hbar = 1.054*10^(-34);%eV*s
% lambda = sqrt(-2*m*E/(hbar^2));
% k = sqrt(2*m*(E-U)/(hbar^2));
% ee = exp(-lambda*L/2);
% ss = sin(k*L/2);
% cc = cos(k*L/2);
% A1 = 1;
% A2 = A1*(1+(tan(k*L/2))^2)/(2*(lambda/k)*tan(k*L/2)+...
%     1-(tan(k*L/2))^2);
% C = (ee*(A1+A2))/(2*cc);
% B = (ee*(A2-A1))/(2*ss);
% xi = 11*10^(-10);
% x = (xi-(L/2):0.0001*10^(-10):xi+(L/2));
% 
% f = B*sin(k*(x-xi))+C*cos(k*(x-xi));
% si = abs(f);
% si = si.^2;
% 
% figure(4)
% plot(x,si)