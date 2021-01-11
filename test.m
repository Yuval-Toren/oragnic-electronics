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
x = (6*10^(-10):0.001*10^(-10):14*10^(-10));

outLeft=exp(lambda.*(x-xi));

fun = @(x) (exp(lambda.*(x-xi))).^2;
N1 = sqrt(integral(fun,0,xi-(L/2)))
a = trapz(x,((1/3.8005e-06)*exp(lambda.*(x-xi))).^2)
outRight=exp(-lambda.*(x-xi));
in=(B.*sin(k.*(x-xi))+C.*cos(k.*(x-xi)));

phi1_1=outLeft(find(x<(xi-L/2)));
phi1_2=in(find(x>=(xi-L/2) & x<=(xi+L/2)));
phi1_3=outRight(find(x>(xi+L/2)));
phi1=[phi1_1,phi1_2,phi1_3];
figure(3)
plot(x,phi1.^2);