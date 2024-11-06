% BOIF1302
% DESM1210

clc
clear all
close all

%% Fonctions de transferts
TF_AZ_num = [1.59e09];
TF_AZ_den = [1 1020.51 25082.705 3102480.725 64155612.5 82700000 0];

TF_EL_num = [7.95e09];
TF_EL_den = [1 1020.51 37082.705 15346520.725 320776412.5 413500000 0];

TF_AZ = tf(TF_AZ_num, TF_AZ_den)
TF_EL = tf(TF_EL_num, TF_EL_den)

%% Variables Télescope A
Mp = 25; %En pourcent
Ts = 1; %sec
Tr_10 = 0.18; %sec

%Erreurs
ERP_unitaire_AZ = 0;
ERP_unitaire_EL = 0;
ERP_rampe_AZ = 0.03; %deg
ERP_rampe_EL = 0; %deg
ERP_para_EL = 0.08; %deg

%% Variables Télescope B


%% Conception Télescope A
Phi = rad2deg(atan(-pi/log(Mp/100)));
Zeta = cosd(Phi);

%On doit trouver le plus grand des Omega_n
% Omega_n = (1+(1*Zeta)+(1.4*(Zeta^2)))/Tr_10;
Omega_n = 4/(Ts*Zeta);

%On trouve Omega_a pour simplifier les P_etoile
Omega_a = Omega_n*sqrt(1-Zeta^2);

%On trouve P_etoile
P_etoile_AZ = -Zeta*Omega_n + Omega_a*i;


%% Calcul pour azimut Télescope 1
frsp = evalfr(TF_AZ, P_etoile_AZ);
Angle_AZ = (rad2deg(angle(frsp)));
clear frsp

%Valider c'est de combien qu'on ajouter
% figure
% hold on
% pzmap(TF_AZ)
% scatter(real(P_etoile), imag(P_etoile), 100, "square", 'black')
%En validant avec le pzmap il est possible de voir qu'on dépassera pas 360 degree

%Calculs qui aideront dans les étapes suivantes
Angle_AZ = Angle_AZ - 360;
Delta_Phi_AZ = - 180 - Angle_AZ;
Alpha_AZ = 180 - Phi;

%Trouver les aangles des distances
Phi_Z_AZ = (Alpha_AZ+Delta_Phi_AZ)/2;
Phi_P_AZ = (Alpha_AZ-Delta_Phi_AZ)/2;

%Trouver les poles et zeros
Z_AZ = real(P_etoile_AZ)-(imag(P_etoile_AZ)/tand(Phi_Z_AZ));
P_AZ = real(P_etoile_AZ)-(imag(P_etoile_AZ)/tand(Phi_P_AZ));

%Cree une sous fonction de transfert pour trouver le Ka
TF_AZ_Ka = tf([1 -Z_AZ], [1 -P_AZ]) * TF_AZ;

%Calcul du Ka
K_AvPh_AZ = 1/abs(evalfr(TF_AZ_Ka, P_etoile_AZ));


%% Variables Télescope B