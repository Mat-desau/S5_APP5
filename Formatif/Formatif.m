 close all
 clear all
 clc

%% Problème 1

% Question 1
%Non cela rends plus stable

%Question 2
%Oui cela peut devenir instable

%Question 3
%Faux car PM = sqrt(tan(PM)sin(PM)

%Question 4
%Faux ca modifie le transitoire

%Question 5
%Faux c'est avec un avance

%Question 6
%Vrai

%Question 7
%Vrai car la bande passante c'est avant que ca traverse 0

%Question 8
Zeta = 0.5*sqrt(tand(60)*sind(60));
%Vrai

%Question 9
%Vrai les deux sont positif

%Question 10
%Vrai car t'es au dessus de -180 degrée

%Question 11
%Faux le pole est a l'infini

%Question 12
%Vrai parce que tu amplifies les hautes fréquences

%Question 13
%Faux ce serait un passe bas

%Question 14
%Vrai car tu montes ta marge de gain donc ton 0 se déplace dans ton bode

%Question 15
%Vrai tu peux prendre des mesures sur ton système sans savoir sa fonction
%de transfert et quand même l'analyser

%% Problème 2
clear all
close all
clc

G = tf([1], [1 2]) * tf([1], [1 2]) * tf([1], [1 3])

%Etape 1 trouver P*
%Etape 2 trouver angle
%Etape 3 Delta phi
%Etape 4 Angle poles zero
%Etape 5 Gain

Phi = atand((-pi)/log(25/100))

Zeta = cosd(Phi)

Wn = 4/(Zeta*1.6)

Wa = Wn*sqrt(1-Zeta^2)

%a) 
P_etoile = (-Zeta*Wn) + (Wa*i)

%b)
Angle = rad2deg(angle(evalfr(G, P_etoile)))-360
Delta_Phi = -180 - Angle

Coordo = P_etoile - (-1)
Phi_z = rad2deg(angle(Coordo))

Variation = Phi_z - Delta_Phi

%c)
Phi_p = -(Delta_Phi - Phi_z)

P = real(P_etoile)-(imag(P_etoile)/tand(Phi_p));
Z = -1

%d)
TF = tf([1 -Z], [1 -P]);
temp = evalfr(TF, P_etoile);
temp2 = evalfr(G, P_etoile);

Ka = 1 /(abs(temp) * abs(temp2));

%e)
TF_Finale = Ka * TF * G;
TF_Finale_BF = feedback(TF_Finale, 1);

% figure
% hold on
% rlocus(TF_Finale, "red")
% plot(real(P_etoile), imag(P_etoile), 'color', "blue", 'marker', "pentagram")
% plot(real(P_etoile), -imag(P_etoile), 'color', "blue", 'marker', "pentagram")

%C'est bon ils sont présent

%f)
% figure
% step(TF_Finale_BF)
stepinfo(TF_Finale_BF);

%Ce ne sont pas bon puisque settlingTime trop haut et Overshoot aussi

%g)
P_etoile = (-Zeta*Wn) + (Wa*i);

Angle = rad2deg(angle(evalfr(G, P_etoile)))-360;
Delta_Phi = (-180 - Angle)/2;

Coordo = P_etoile - (-1);
Phi_z = rad2deg(angle(Coordo));

Variation = Phi_z - Delta_Phi;

Phi_p = -(Delta_Phi - Phi_z);

P = real(P_etoile)-(imag(P_etoile)/tand(Phi_p));
Z = -1;

TF = tf([1 -Z], [1 -P]) * tf([1 -Z], [1 -P]);
temp = evalfr(TF, P_etoile);
temp2 = evalfr(G, P_etoile);

Ka = 1 /(abs(temp) * abs(temp2));

TF_Finale = Ka * TF * G
TF_Finale_BF = feedback(TF_Finale, 1);

% figure
% hold on
% rlocus(TF_Finale, "red")
% plot(real(P_etoile), imag(P_etoile), 'color', "blue", 'marker', "pentagram")
% plot(real(P_etoile), -imag(P_etoile), 'color', "blue", 'marker', "pentagram")

% figure
% step(TF_Finale_BF)
stepinfo(TF_Finale_BF);

%% Probleme 3
clear all 
close all
clc

%a)
G = tf([1], [1 5 0]) * tf([1], [1 11]);

Kp = 218.7;

TF = Kp * G;
TF_BF = feedback(TF, 1);

[num, den] = tfdata(TF_BF, 'v');
[num2, den2] = tfdata(TF, 'v');

[r, p, k] = residue(num, den);

poids = abs(r)./abs(real(p));

TF2 = tf([r(2)], [1 -p(2)]) + tf([r(3)], [1 -p(3)]);

Poles = roots(den)

%ts=4/zeta*wn et on sait que le real du pole c'est wn
Ts = 4/(-real(Poles(2)))

Phi = 180 - rad2deg(angle(Poles(2)));
Mp = 100*exp(-pi/tand(Phi))

K = num2(end)/den2(end-1)
Erreur = 1 / K

%b)


%% Problème 4 
clear all
close all
clc

G = tf([1], [1 8 0]) * tf([1], [1 30]);

Phi = atand((-pi)/log(10/100));

Zeta = cosd(Phi);

Wa = pi/0.4;

Wn = Wa/(sqrt(1-Zeta^2));

%On trouve les désirer
Pm_des = atand((2*Zeta)/(sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2)));
Wg_des = Wn*sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2);

%Calcul du K_etoile
K_etoile = 1 / abs(evalfr(G, Wg_des*i));

% figure
% margin(K_etoile*G)

PM = rad2deg(angle(evalfr(K_etoile*G, Wg_des*i)))+180;
Delta_phi = Pm_des - PM + 5;

Alpha = (1 - sind(Delta_phi))/(1 + sind(Delta_phi));

T = 1/(Wg_des*sqrt(Alpha));

Ka = K_etoile/sqrt(Alpha);

P = -1/T;
Z = -1/(T*Alpha);

TF = tf([1 -P], [1 -Z]);
TF = Ka * TF;
TF_test = TF*G;
[num, den] = tfdata(TF_test, 'v');

% figure 
% margin(TF)

Kvel = num(end)/den(end-1);
Kvel_etoile = 1/0.05;

K_etoile = Kvel_etoile/Kvel;

T = 10/Wg_des;

Z = -1/T;
P = -1/(T*K_etoile);

TF2 = tf([1 -Z], [1 -P]);

TF_Finale = TF * TF2 * G;

% figure
% margin(TF_Finale)

%% Problème 5

clear all
close all
clc

G = tf([1 2 1], [1 2 2 1 0])

ERP = 0.005
PM_des = 60
BW = 20

%Trouver les désirer
Zeta = (0.5) * sqrt(tand(PM_des)*sind(PM_des))
Wg_des = BW * sqrt(sqrt(1+4*Zeta^4)-(2*Zeta^2))/(sqrt((1-2*Zeta^2)+sqrt((4*Zeta^4)-(4*Zeta^2)+2)))

%Trouver K*
K_etoile = 1/abs(evalfr(G,Wg_des*i))

% figure
% margin(G*K_etoile)

%Compensation PM
Angle = rad2deg(angle(evalfr(G, Wg_des*i)))
PM = Angle + 180

%On trouve delta Phi
Delta_phi = PM_des - PM + 5

%Trouver Alpha
Alpha = (1 - sind(Delta_phi))/(1 + sind(Delta_phi))

% On trouves poles et zeros
T = 1/(Wg_des*sqrt(Alpha))

Z = -1/T
P = -1/(T*Alpha)

%Trouver Ka
Ka = K_etoile / sqrt(Alpha)

TF = Ka*tf([1 -Z], [1 -P])
TF_Finale_Av = TF * G
[num, den] = tfdata(TF_Finale_Av, 'v')

Kvel = num(end)/den(end-1)
Kvel_etoile = 1/ERP

%On trouve Ketoile
K_etoile = Kvel_etoile / Kvel

%On trouve les poles et zeros
T = 10/Wg_des

Z = -1/T
P = -1/(K_etoile*T)

TF = tf([1 -Z], [1 -P])
TF_Finale = Ka*TF*TF_Finale_Av

%% Problème 6
clear all
close all
clc

G = tf([1 2 1], [1 2 2 1 0]);

Phi = atand(-pi/log(10/100));

Zeta = cosd(Phi);

Wn = 4/(1*Zeta);

ERP = 0.005;

Pm_des = atand((2*Zeta)/(sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2)));
Wg_des = Wn*sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2);

K_etoile = 1/abs(evalfr(G, Wg_des*i));

% figure
% margin(K_etoile*G)

Angle = rad2deg(angle(evalfr(K_etoile*G, Wg_des*i)));
PM = Angle + 180;

Delta_phi = Pm_des - PM;

Alpha = (1 - sind(Delta_phi))/(1 + sind(Delta_phi));

T = 1/(Wg_des*sqrt(Alpha));

Z = -1/T;
P = -1/(Alpha*T);

Ka = K_etoile / sqrt(Alpha);

TF = tf([1 -Z], [1 -P]);
TF_Finale_Av = Ka * TF * G;
[num, den] = tfdata(TF_Finale_Av, 'v');

Kvel = num(end)/den(end-1);
Kvel_etoile = 1/ERP;

K_etoile = Kvel_etoile/Kvel;

T = 10/Wg_des;

Z = -1/T;
P = -1/(K_etoile*T);

TF = tf([1 -Z], [1 -P]);

TF_Finale = TF_Finale_Av * TF;
[num, den] = tfdata(TF_Finale, 'v');

K_trouver = num(end)/den(end-1);

EPS_trouver = 1 / K_trouver;

%% Probleme 7
clear all
close all
clc

G = tf([1], [1 8 0]) * tf([1], [1 30])

K = 1075

ERP = 0.5

TF_Finale_K = K * G
TF_Finale_BF = feedback(TF_Finale_K, 1)

[num, den] = tfdata(TF_Finale_BF, 'v');
[num2, den2] = tfdata(TF_Finale_K, 'v');

Kacc_des = 1 / ERP
Kvel_now = num2(end)/den2(end-1)

%a)
% figure
% margin(TF_Finale)

[~,~,~,Wn] = margin(TF_Finale_K)
Poles = roots(den)

%b)
Z = real(Poles(2))/10

Ki = Kacc_des / Kvel_now
Kp = -Ki/Z

TF = Kp*tf([1 -Z], [1 0])

TF_Finale_PI = TF * TF_Finale_K

%c)
Wn = abs(Poles(2))
Zeta = real(Poles(2))/Wn

Wg_des = Wn*sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2)

Z = -Wg_des/10

Ki = Kacc_des / Kvel_now

Kp = -Ki/Z

TF = Kp*tf([1 -Z], [1 0])

TF_Finale_PI = TF * TF_Finale_K