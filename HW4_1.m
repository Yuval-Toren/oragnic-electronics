%% 1 a 1
clear;clc;
U = -1.922*(10^(-18));%eV
E = (U:0.001*(10^(-18)):0);%eV
m = 9.10*10^(-31);%Kg
L = 1*10^(-10);%m
hbar = 1.054*10^(-34);%eV*s
lambda = sqrt(-2*m*E/(hbar^2));
k = sqrt(2*m*(E-U)/(hbar^2));
%tan(sqrt(k)*L) = (2*sqrt(lambdas.*ks))/(ks.^2-lambdas.^2);
figure(1)
plot(E,tan(k*L))
hold on
plot(E,(2*lambda.*k)./(k.^2-lambda.^2))
xlabel('E[J]');
ylabel('[a.u.]');
title('-12eV Well');
legend('tan(kL)','(2lambda*k)/(k^2-lambda^2)')
%% 1 a 2
U = -2.035*(10^(-18));%eV
E = (U:0.001*(10^(-18)):0);%eV
m = 9.10*10^(-31);%Kg
L = 1.1*10^(-10);%m
hbar = 1.054*10^(-34);%eV*s
lambda = sqrt(-2*m*E/(hbar^2));
k = sqrt(2*m*(E-U)/(hbar^2));

figure(2)
plot(E,tan(k*L))
hold on
plot(E,(2*lambda.*k)./(k.^2-lambda.^2))
xlabel('E[J]');
ylabel('[a.u.]');
title('-12.7eV Well');
legend('tan(kL)','(2lambda*k)/(k^2-lambda^2)')
%% 1 b 1
U = -1.922*(10^(-18));%eV
E = -7.718*(10^(-19));%eV
m = 9.10*10^(-31);%Kg
L = 1*10^(-10);%m
hbar = 1.054*10^(-34);%eV*s
lambda = sqrt(-2*m*E/(hbar^2));
k = sqrt(2*m*(E-U)/(hbar^2));
ee = exp(-lambda*L/2);
ss = sin(k*L/2);
cc = cos(k*L/2);
mat = [ee,0,ss,-cc;...
    0,ee,-ss,-cc;...
    lambda*ee*10^(-10),0,-k*cc*10^(-10),-k*ss*10^(-10);...
    0,-lambda*ee*10^(-10),-k*cc*10^(-10),k*ss*10^(-10);];
A1 = 1;
A2 = A1*(1+(tan(k*L/2))^2)/(2*(lambda/k)*tan(k*L/2)+...
    1-(tan(k*L/2))^2);
C = (ee*(A1+A2))/(2*cc);
B = (ee*(A2-A1))/(2*ss);
xi = 8*10^(-10);
x = (0:0.001*10^(-10):16*10^(-10));

outLeft=A1.*exp(lambda.*(x-xi));
outRight=A2*exp(-lambda.*(x-xi));
in=(B.*sin(k.*(x-xi))+C.*cos(k.*(x-xi)));

phi1_1=outLeft((x<(xi-L/2)));
phi1_2=in((x>=(xi-L/2) & x<=(xi+L/2)));
phi1_3=outRight((x>(xi+L/2)));
phi1=[phi1_1,phi1_2,phi1_3];
Nphi1 = phi1/sqrt(trapz((x(2)-x(1))*abs(phi1).^2));
figure(3)
plot(x,abs(Nphi1).^2);
xlabel('x[m]');
ylabel('Probability Density');
title('-12eV Well');
%% 1 b 2
U = -2.035*(10^(-18));%eV
E = -9.1324e-19;%eV
m = 9.10*10^(-31);%Kg
L = 1.1*10^(-10);%m
hbar = 1.054*10^(-34);%eV*s
lambda = sqrt(-2*m*E/(hbar^2));
k = sqrt(2*m*(E-U)/(hbar^2));
ee = exp(-lambda*L/2);
ss = sin(k*L/2);
cc = cos(k*L/2);
mat = [ee,0,ss,-cc;...
    0,ee,-ss,-cc;...
    lambda*ee*10^(-10),0,-k*cc*10^(-10),-k*ss*10^(-10);...
    0,-lambda*ee*10^(-10),-k*cc*10^(-10),k*ss*10^(-10);];
A1 = 1;
A2 = A1*(1+(tan(k*L/2))^2)/(2*(lambda/k)*tan(k*L/2)+...
    1-(tan(k*L/2))^2);
C = (ee*(A1+A2))/(2*cc);
B = (ee*(A2-A1))/(2*ss);
xi = 11*10^(-10);
x = (0:0.001*10^(-10):16*10^(-10));

outLeft=A1.*exp(lambda.*(x-xi));
outRight=A2*exp(-lambda.*(x-xi));
in=(B.*sin(k.*(x-xi))+C.*cos(k.*(x-xi)));

phi2_1=outLeft((x<(xi-L/2)));
phi2_2=in((x>=(xi-L/2) & x<=(xi+L/2)));
phi2_3=outRight((x>(xi+L/2)));
phi2=[phi2_1,phi2_2,phi2_3];
Nphi2 = phi2/sqrt(trapz((x(2)-x(1))*abs(phi1).^2));
figure(4)
plot(x,abs(Nphi2).^2);
xlabel('x[m]');
ylabel('Probability Density');
title('-12.7eV Well');
%% part b

U = -1.922*(10^(-18));%eV
E = -7.718*(10^(-19));%eV
m = 9.10*10^(-31);%Kg
L = 1*10^(-10);%m
hbar = 1.054*10^(-34);%eV*s
lambda = sqrt(-2*m*E/(hbar^2));
k = sqrt(2*m*(E-U)/(hbar^2));
ee = exp(-lambda*L/2);
ss = sin(k*L/2);
cc = cos(k*L/2);
mat = [ee,0,ss,-cc;...
    0,ee,-ss,-cc;...
    lambda*ee*10^(-10),0,-k*cc*10^(-10),-k*ss*10^(-10);...
    0,-lambda*ee*10^(-10),-k*cc*10^(-10),k*ss*10^(-10);];
A1 = 1;
A2 = A1*(1+(tan(k*L/2))^2)/(2*(lambda/k)*tan(k*L/2)+...
    1-(tan(k*L/2))^2);
C = (ee*(A1+A2))/(2*cc);
B = (ee*(A2-A1))/(2*ss);
xi = 11*10^(-10);
x = (0:0.001*10^(-10):16*10^(-10));

outLeft=A1.*exp(lambda.*(x-xi));
outRight=A2*exp(-lambda.*(x-xi));
in=(B.*sin(k.*(x-xi))+C.*cos(k.*(x-xi)));

phi3_1=outLeft((x<(xi-L/2)));
phi3_2=in((x>=(xi-L/2) & x<=(xi+L/2)));
phi3_3=outRight((x>(xi+L/2)));
phi3=[phi3_1,phi3_2,phi3_3];
Nphi3 = phi3/sqrt(trapz((x(2)-x(1))*abs(phi1).^2));


U1=-12; 
U2=-12; 
S11=1;
S22=1;
f=phi1.*phi3;
S12=trapz(x,f);
S21=S12;

E1=-4.8; 
E2=-4.8; 

S=[S11 S22 S12 S21];
E=[E1 E2 E2 E1];
U=[U1 U2 U1 U2];
H=S.*(E-U);
a1 = 1;
a2 = a1;
Phi_sum=a1*phi1+a2*phi3;

figure(5)
plot(x,Phi_sum);
xlabel('x (m)');
ylabel('Wave Function');
title('Identical Wells');


figure(6)
plot(x,abs(Phi_sum).^2);
xlabel('x (m)');
ylabel('Probability Density');
title('Identical Wells');

%%
U1=-12; 
U2=-12.7; 
S11=1;
S22=1;
f=phi1.*phi2;
S12=trapz(x,f);
S21=S12;

E1=-4.8; 
E2=-5.7; 

S=[S11 S22 S12 S21];
E=[E1 E2 E2 E1];
U=[U1 U2 U1 U2];
H=S.*(E-U);
a1 = 1;
a2 = -1/7;
Phi_sum=a1*phi1+a2*phi2;

figure(7)
plot(x,Phi_sum);
xlabel('x (m)');
ylabel('Wave Function');
title('Different Wells');


figure(8)
plot(x,abs(Phi_sum).^2);
xlabel('x (m)');
ylabel('Probability Density');
title('Different Wells');

