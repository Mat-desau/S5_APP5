%Laboratoire 2
close all
clear all
clc
%% Problème 8

G = tf([1 4], [1 8 3 0]);
PM_etoile = 50;
BW_etoile = 5.2;
ERP_Rampe = 0.005;

Zeta = 0.5*sqrt(tand(PM_etoile)*sind(PM_etoile));

Omega_g = BW_etoile * (sqrt(sqrt(1+(4*Zeta^4))-(2*Zeta^2))/(sqrt((1-(2*Zeta^2))+sqrt((4*Zeta^4)-(4*Zeta^2)+2))));

K_etoile = 1 / abs(evalfr(G, (Omega_g*i)));

PM = rad2deg(angle(evalfr((K_etoile*G), (Omega_g*i)))) + 180;

Delta_Phi = PM_etoile - PM + 5;

Alpha = (1 - sind(Delta_Phi))/(1+sind(Delta_Phi));

T = 1 / (Omega_g * sqrt(Alpha));

Z = -1/T;
P = -1/(Alpha*T);

K_a = K_etoile/sqrt(Alpha);

Ga = K_a * tf([1 -Z], [1 -P]);
G_Av = Ga * G;

G_Av_BF = feedback(G_Av, 1);

[num, den] = tfdata(G_Av, 'v');

BW_trouver = bandwidth(G_Av_BF);

Kvel_etoile = 1 / ERP_Rampe;
Kvel = num(end)/den(end-1);

K_etoile = Kvel_etoile/Kvel;

T = 10/(Omega_g);

Z = -Omega_g/10;
P = -1/(K_etoile*T);

Gr = tf([1 -Z], [1 -P]);
G_Av_Re = Gr * G_Av;
G_Av_Re_BF = feedback(G_Av_Re, 1);

BW_trouver = bandwidth(G_Av_Re_BF);

% figure
% margin(G_Av_Re)

%% Problème 9
close all
clear all
clc

G = tf([1], [1 0]) * tf([1], [1 1]) * tf([1], [1 4]);

Mp = 12;
Phi = atand(-pi/log(Mp/100));

Zeta = cosd(Phi);

PM_etoile = atand((2*Zeta)/sqrt(sqrt(1+(4*Zeta^4))-(2*Zeta^2)));

phase = (-1) * (180-PM_etoile);

% figure
% margin(G)

Omega_g = 0.499;

K = 1 / abs(evalfr(G, Omega_g*i));

% [mag, pha] = bode(G, Omega_g)
% K_etoile = 1/mag;

KG = K * G;

[num, den] = tfdata(KG, 'v');

Kvel = num(end)/den(end-1);





disp("Hello World")