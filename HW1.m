%% Homework 1 - Yuval Toren 303031603
%% Problem 3
%% J(V)
J = (0:0.01*10:1*10);
mu = 10^-9;
L = 75*10^-9;
x = (0:0.01*10^-9:L);
n = 1.8;
epsr = n^2;
N0 = 5*10^21;
eps0 = 8.85*10^-12;
q = 1.6*10^-19;%1.6*10^-19;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
V = sqrt((8.*J)/(9*epsr*eps0*mu)).*((K+L).^(3/2)-K.^(3/2));
%% plots
f1 = figure(1);
plot(V,J);
title("J VS V");
xlabel("V[V]")
ylabel("J[A/m^2]");
hold on
N0 = 1*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
V = sqrt((8.*J)/(9*epsr*eps0*mu)).*((K+L).^(3/2)-K.^(3/2));
plot(V,J);
hold on
N0 = 5*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
V = sqrt((8.*J)/(9*epsr*eps0*mu)).*((K+L).^(3/2)-K.^(3/2));
plot(V,J);
legend('N_0 = 5*10^2^1 [m^-^3]','N_0 = 1*10^2^1 [m^-^3]','N_0 = 5*10^2^2 [m^-^3]');
saveas(f1,"JV.jpg")
%% E(X)
N0 = 5*10^21;
J = 5;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
E = (2*J/(epsr*eps0*mu))*(x+K);
%% plots
f2 = figure(2);
plot(x,E);
title("E VS x");
xlabel("x[m]")
ylabel("E[V/m]");
hold on
N0 = 1*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
E = (2*J/(epsr*eps0*mu))*(x+K);
plot(x,E);
hold on
N0 = 5*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
E = (2*J/(epsr*eps0*mu))*(x+K);
plot(x,E);
legend('N_0 = 5*10^2^1 [m^-^3]','N_0 = 1*10^2^1 [m^-^3]','N_0 = 5*10^2^2 [m^-^3]');
saveas(f2,"Ex.jpg")
%% n(x)
N0 = 5*10^21;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
N = sqrt(J*epsr*eps0/(2*mu))*(1/q)*(1./sqrt(K+x));
%% plots
f5 = figure(5);
plot(x,N);
title("n VS x");
xlabel("x[m]")
ylabel("n[1/m^3]");
hold on
N0 = 1*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
N = sqrt(J*epsr*eps0/(2*mu))*(1/q)*(1./sqrt(K+x));
plot(x,N);
hold on
N0 = 5*10^22;
K = (epsr*eps0*J)/(2*mu*q^2*N0^2);
N = sqrt(J*epsr*eps0/(2*mu))*(1/q)*(1./sqrt(K+x));
plot(x,N);
legend('N_0 = 5*10^2^1 [m^-^3]','N_0 = 1*10^2^1 [m^-^3]','N_0 = 5*10^2^2 [m^-^3]');
saveas(f5,"Nx.jpg")
%% Problem 4
%% V(x)
epsr = 3^2;
L = 100*10^(-9);
Nc = 10^27;
KT = 298*8.62*10^(-5);
x = (0:0.01*10^-9:L);
V = (q*Nc/(epsr*eps0))*((L*KT)/0.4)*exp(-0.3/KT).*...
    ((1+exp((-0.4.*x)./(L*KT))).*((L*KT)/0.4)+x);
%% plots
f4 = figure(4);
plot(x,V);
title("V VS x");
xlabel("x[m]")
ylabel("V[V]");
saveas(f4,"Vx.jpg")