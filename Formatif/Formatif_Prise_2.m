%% Problème 2
clc
close all
clear all

G = tf([1], [1 2]) * tf([1], [1 2]) * tf([1], [1 3])

Ts = 1.6
Mp = 25

% Trouver les valeurs
Phi = atand(-pi/(log(Mp/100)))

Zeta = cosd(Phi)

Wn = 4 / (Zeta*Ts)

Wa = Wn*sqrt(1-Zeta^2)

%Trouver P_etoile
P_etoile = (-Zeta*Wn) + Wa*i

Z = -1

Angle = rad2deg(angle(evalfr(G, P_etoile))) - 360
Delta_phi = -180 - Angle

Coordo = P_etoile - Z
Phi_z = rad2deg(angle(Coordo))

Phi_p = -Delta_phi + Phi_z

P = real(P_etoile)-(imag(P_etoile)/tand(Phi_p))

TF = tf([1 -Z], [1 -P])
temp1 = evalfr(TF, P_etoile)
temp2 = evalfr(G, P_etoile)

Ka = 1 / abs(temp1*temp2)

FT_Comp = TF*Ka

FT_Finale = FT_Comp*G
FT_Finale_BF = feedback(FT_Finale, 1)
% figure
% rlocus(FT_Finale, "red")
stepinfo(FT_Finale_BF)

%% Probleme 4
clc
close all
clear all

G = tf([1], [1 8 0]) * tf([1], [1 30])

Tp = 0.4
Mp = 10
ERP = 0.05

Phi = atand(-pi/(log(Mp/100)))
Zeta = cosd(Phi)

Wa = pi / Tp

Wn = Wa / sqrt(1-Zeta^2)

Pm_des = atand((2*Zeta)/(sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2)));
Wg_des = Wn*sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2);

K_etoile = 1 / abs(evalfr(G, Wg_des*i))

% figure 
% margin(K_etoile*G)

PM = rad2deg(angle(evalfr(K_etoile*G, Wg_des*i))) + 180 

Delta_phi = Pm_des-PM

Alpha = (1-sind(Delta_phi))/(1+sind(Delta_phi))

T = 1/(Wg_des*sqrt(Alpha))

Z = -1/T
P = -1/(Alpha*T)

Ka = K_etoile/sqrt(Alpha)

TF = tf([1 -Z], [1 -P])

TF_comp = TF*Ka
TF_Finale = TF_comp * G
[num, den] = tfdata(TF_Finale, 'v')

% Retard

Kvel_Cal = num(end)/den(end-1)
Kvel_des = 1 / ERP

K_etoile = Kvel_des/Kvel_Cal

T = 10/Wg_des

Z = -1/T
P = -1/(K_etoile*T)

TF = tf([1 -Z], [1 -P])

TF_Finale = TF * TF_Finale
TF_Finale_BF = feedback(TF_Finale, 1)
[num, den] = tfdata(TF_Finale, 'v')

Kvel = 1 / (num(end)/den(end-1))

stepinfo(TF_Finale_BF)
