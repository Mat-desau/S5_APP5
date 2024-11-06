clear all
close all
clc

%% Problème 5
Mp = 10; % en %
ts = 4;

num = [2];
den = [2 5 2]; 

G = tf(num, den) * tf([1], [1 0]);

Phi = atan(-pi/log(Mp/100));

Zeta = cos(Phi);

Omega_n = 4/(ts*Zeta);

Omega_a = Omega_n*sqrt(1-Zeta^2);

P_etoile = -Zeta*Omega_n + Omega_a*i;

frsp = evalfr(G, P_etoile);
Angle = (rad2deg(angle(frsp))-360);

Delta_Phi = - 180 - Angle;

Z = real(P_etoile)-(imag(P_etoile)/tand(Delta_Phi));

TF1 = tf([1 -Z], [1]);
TF2 = TF1*G;

Kp = 1/abs(evalfr(TF2, P_etoile));

% figure
% hold on
% rlocus(num, den, "red")
% scatter(Z, 0, 50, "blue")

%b)
G2 = TF2*Kp;

Kvel = 2.862/2;

E = 1/Kvel;

%c)
Mp = 10; % en %
ts = 4;

num = [2];
den = [2 5 2]; 

G = tf(num, den) * tf([1], [1 0]);

Phi = atan(-pi/log(Mp/100));

Zeta = cos(Phi);

Omega_n = 4/(ts*Zeta);

Omega_a = Omega_n*sqrt(1-Zeta^2);

P_etoile = -Zeta*Omega_n + Omega_a*i;

frsp = evalfr(G, P_etoile);
Angle = (rad2deg(angle(frsp))-360);

Delta_Phi = (- 180 - Angle)/2;

Z = real(P_etoile)-(imag(P_etoile)/tand(Delta_Phi));

TF1 = tf([1 -Z], [1]);
TF2 = (TF1^2)*G;

Kp = 1/abs(evalfr(TF2, P_etoile));

G2 = TF2*Kp;

Kvel = 11.45/2;

E = 1/Kvel;

%% Probleme 6
clear all
close all
clc

Mp = 6;
Tr10 = 0.004;
Tp = 0.008;
Ts = 0.010;
erp = 0.00005;

G = tf([4500], [1 361.2 0]);

%a)
Phi = atan(-pi/log(Mp/100));
Zeta = cos(Phi);

%Omega_n = (1+(1*Zeta)+(1.4*(Zeta^2)))/Tr10
Omega_n = 4/(Ts*Zeta);
%Omega_n = pi/ (Tp * (sqrt(1-(Zeta^2))))

Omega_a = Omega_n*sqrt(1-Zeta^2);

P_etoile = -Zeta*Omega_n + Omega_a*i;

%b)
frsp = evalfr(G, P_etoile);
Angle = rad2deg(angle(frsp))-360;

Delta_Phi = -180 - Angle;
Alpha = 180 - rad2deg(Phi);

Phi_z = (Alpha + Delta_Phi)/2;
Phi_p = (Alpha - Delta_Phi)/2;

Za = real(P_etoile) - imag(P_etoile)/tand(Phi_z);
Pa = real(P_etoile) - imag(P_etoile)/tand(Phi_p);

G_a = tf([1 -Za], [1 -Pa]) * G;

Ka = 1/abs(evalfr(G_a, P_etoile));

G_a = G_a * Ka;

%c)
Kvel = (1.612e08)/(3.204e05);
Kvel_etoile = 1/erp;
K_etoile = Kvel_etoile/Kvel;

Zr = real(P_etoile)/10;
Pr = Zr/K_etoile;

G_r = tf([1 -Zr], [1 -Pr]) * G_a;

Kr = 1/abs(evalfr(G_r, P_etoile));

G_r = G_r * Kr;

%e)
Zpd = real(P_etoile) - imag(P_etoile)/tand(Delta_Phi);

G_pd = tf([1 -Zpd], [1]) * G;

Kpd = 1/abs(evalfr(G_pd, P_etoile));

G_pd = G_pd * Kpd;

%f)
Zpi = real(P_etoile)/10;

G_pi = G_pd *tf([1 -Zpi], [1 0]);

%% Problème 7
clc
clear all
close all

G = tf([1 1], [1 6 10]) * tf([1], [1 7]) * tf([1], [1 0])

Ts = 0.1;
Mp = 5;

Phi = rad2deg(atan(-pi/log(Mp/100)))
Zeta = cosd(Phi)

Omega_n = 4/(Ts*Zeta)
Omega_a = Omega_n*sqrt(1-Zeta^2)

P_etoile = -Zeta*Omega_n + Omega_a*i

frsp = evalfr(G, P_etoile)
Angle = rad2deg(angle(frsp)) - 360

Delta_Phi = -180 - Angle
Alpha = 180 - Phi

Zpd = real(P_etoile) - imag(P_etoile)/tand(Delta_Phi/2)

G_pd2 = tf([1 -Zpd], [1]) * tf([1 -Zpd], [1]) * G
