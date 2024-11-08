% BOIF1302
% DESM1210
clc
clear all
close all

%% Ajout de la trajectoire
Profile_Tracking
%la matrice-colonne de temps en secondes est dans ttrk
%la matrice-colonne de consigne en degrés est dans utrk

%% Fonctions de transferts
TF_AZ_num = [1.59e09];
TF_AZ_den = [1 1020.51 25082.705 3102480.725 64155612.5 82700000 0];

TF_EL_num = [7.95e09];
TF_EL_den = [1 1020.51 37082.705 15346520.725 320776412.5 413500000 0];

TF_AZ = tf(TF_AZ_num, TF_AZ_den);
TF_EL = tf(TF_EL_num, TF_EL_den);

%% Variables Télescope A
Mp_A = 25; %En pourcent
Ts_A = 1; %sec
Tr_10_A = 0.18; %sec

%Erreurs
ERP_unitaire_A = 0;
ERP_rampe_AZ_A = 0.03; %deg
ERP_rampe_EL_A = 0; %deg
ERP_para_EL_A = 0.08; %deg

%% Variables Télescope B
BW_B = 10; %rad/s
PM_B = 50; %deg +- 1 deg

%Erreurs
ERP_rampe_B = 0.005; %deg
t_ERP_rampe_B = 8; %sec

%% Conception spécifications Télescope A

disp("_______________________________________ AZ ________________________________")
Phi = rad2deg(atan(-pi/log(Mp_A/100)));
Zeta = cosd(Phi);

%On doit trouver le plus grand des Omega_n
%Omega_n = (1+(1*Zeta)+(1.4*(Zeta^2)))/Tr_10_A
Omega_n = 4/(Ts_A*Zeta);

%On trouve Omega_a pour simplifier les P_etoile
Omega_a = Omega_n*sqrt(1-Zeta^2);

%Ajustement de P_etoile
Ajout_Omega_n = -1.3;
Ajout_Omega_a = -5;

%On trouve P_etoile
P_etoile_A = (-Zeta*Omega_n + Ajout_Omega_n) + (Omega_a + Ajout_Omega_a)*i;

%% Calcul pour Avance phase Azimut Télescope A
frsp = evalfr(TF_AZ, P_etoile_A);
Angle_AZ = (rad2deg(angle(frsp)));
clear frsp

%Calculs qui aideront dans les étapes suivantes
Angle_AZ = Angle_AZ - 360;
Delta_Phi_AZ = - 180 - Angle_AZ;
Alpha_AZ = 180 - Phi;

%Trouver les angles des distances
Phi_Z_AZ = (Alpha_AZ+Delta_Phi_AZ)/2;
Phi_P_AZ = (Alpha_AZ-Delta_Phi_AZ)/2;

%Trouver les poles et zeros
Z_AZ = real(P_etoile_A)-(imag(P_etoile_A)/tand(Phi_Z_AZ));
P_AZ = real(P_etoile_A)-(imag(P_etoile_A)/tand(Phi_P_AZ));

%Cree une sous fonction de transfert pour trouver le Ka
TF_Ka_AZ = tf([1 -Z_AZ], [1 -P_AZ]) * TF_AZ;

%Calcul du K_AvPh_AZ
K_AvPh_AZ = 1/abs(evalfr(TF_Ka_AZ, P_etoile_A));

%Nouvelle fonction de transfert d'avance de phase
TF_AvPh_AZ = TF_Ka_AZ * K_AvPh_AZ;
 
%% Calcul pour retard phase cascades Azimut Télescope A
Diviser = 10;

%Trouver les valeurs des numérateurs et denominateur
[num_temp, den_temp] = tfdata(TF_AvPh_AZ, 'v');

%Trouver les K_etoile avec erreurs
Kvel_AZ = (num_temp(end))/(den_temp(end-1));
Kvel_etoile_AZ = 1/ERP_rampe_AZ_A;
K_etoile_AZ = Kvel_etoile_AZ/Kvel_AZ;
clear num_temp den_temp

%Trouver poles et zeros
Z_RePh_AZ = real(P_etoile_A)/Diviser;
P_RePh_AZ = Z_RePh_AZ/K_etoile_AZ;

%Cree une sous fonction de transfert pour trouver le Ka
TF_Kr_AZ = tf([1 -Z_RePh_AZ], [1 -P_RePh_AZ]) * TF_AvPh_AZ;

%Calcul du K_RePh_AZ
K_RePh_AZ = 1/abs(evalfr(TF_Kr_AZ, P_etoile_A));
%on voit que c'est environ 1 donc on change pour a
K_RePh_AZ = 1;

%Nouvelle fonction de transfert d'avance de phase et retard
TF_RePh_AvPh_AZ = TF_Kr_AZ * K_RePh_AZ;

%% Coupe bande AZ Téléscope 
Omega_o = 54.8; %Peak sur le bode
X = 0.2; %
Kfcp = 1;

num = Kfcp*[1 1 Omega_o.^2];
den = [1 2*X*Omega_o Omega_o.^2];

TF_Coupe_Bande_AZ = tf(num, den);

clear num den Omega_o X Kfcp

%Mettre Coupe-Bande sur les fonctions de transferts
TF_Finale_AZ = TF_Coupe_Bande_AZ * TF_RePh_AvPh_AZ;
TF_Finale_BF_AZ = feedback(TF_Finale_AZ, 1);


%% Demande pour rapport
%Rlocus du système avec les P desirer
% figure
% hold on
% rlocus(TF_AZ, "red")
% rlocus(TF_Finale_AZ, "blue")
% scatter(real(P_etoile_A), imag(P_etoile_A), 50, "*", 'black')
% scatter(real(P_etoile_A), imag(-P_etoile_A), 50, "*", 'black')
% legend(["Original", "Finale", "Pole1", "Pole2"])
% xlim([-30, 10])
% ylim([-100, 100])

%Réponse à l'échelon unitaire
% figure;
% step(TF_Finale_BF_AZ);

%Erreur à une rampe uniaire du système
temps = [0:0.1:30];
Rampe = [0:0.1:30];
y = lsim(TF_Finale_BF_AZ, Rampe, temps);

% figure
% hold on
% box on
% plot(temps, y', "blue")
% plot(temps, Rampe-y', "red")
% plot(temps, Rampe, "black")
% legend(["Réponse", "Erreur", "Rampe"]);
% title("Erreur à une rampe unitaire");
% ylabel("Amplitude");
% xlabel("Time (secondes)")

%Diagramme de bode du système
% figure
% margin(TF_Finale_AZ)

%Erreur sur la trajectoire
% figure
% hold on
% box on
% y = lsim(TF_Finale_BF_AZ, utrk, ttrk);
% plot(ttrk, y, "blue")
% plot(ttrk, utrk-y, "red")
% plot(ttrk, utrk, "black")
% legend(["Réponse", "Erreur", "Trajectoire"]);
% %legend(["Erreur"])
% title("Erreur sur la trajectoire");
% ylabel("Amplitude");
% xlabel("Time (secondes)")


%% Validation système AZ
%On vérifie Mp < 30%    Tr<0.25sec      ts< 1.20sec
stepinfo(TF_Finale_BF_AZ, RiseTimeThreshold=[0 1]);

%on vérifie GM > 10 dB      RM > 0.09s
[Gm, Pm, wcg, wcp] = margin(TF_Finale_AZ);
Pm;
Gm = 20*log10(Gm);
Rm = (Pm/wcp)*(pi/180);

Z_AZ;
P_AZ;

K_AvPh_AZ;
TF_AvPh_AZ;

Kvel_AZ;
K_etoile_AZ;

P_RePh_AZ;
Z_RePh_AZ;

K_RePh_AZ;
TF_RePh_AvPh_AZ;

%disp("Temps Erreur Rampe");
lsiminfo(Rampe-y', temps);

                    %Effacer les non utiliser
                    clear Gm Pm wcg wcp
                    %Clear les variables
                    clear Kvel_AZ Kvel_etoile_AZ K_etoile_AZ Z_RePh_AZ P_RePh_AZ TF_Kr_AZ
                    clear P_etoile_A Ajout_Omega_n Ajout_Omega_a
                    %Clear les variables
                    clear TF_Ka_AZ P_AZ Z_AZ Phi_P_AZ Phi_Z_AZ Angle_AZ Delta_Phi_AZ Alpha_AZ







disp("_______________________________________ EL ________________________________")
%% Calcul pour avance phase cascades Elevation Télescope A
%Ajustement de P_etoile
Ajout_Omega_n = 0;
Ajout_Omega_a = 0;

%On trouve P_etoile
P_etoile_A = -Zeta*(Omega_n + Ajout_Omega_n) + (Omega_a + Ajout_Omega_a)*i;

frsp = evalfr(TF_EL, P_etoile_A);
Angle_EL = (rad2deg(angle(frsp)));
clear frsp

%Calculs qui aideront dans les étapes suivantes
Angle_EL = Angle_EL - 360;
Delta_Phi_EL = - 180 - Angle_EL;
Alpha_EL = 180 - Phi;

%Trouver les angles des distances
Phi_Z_EL = (Alpha_EL+Delta_Phi_EL)/2;
Phi_P_EL = (Alpha_EL-Delta_Phi_EL)/2;

%Trouver les poles et zeros
Z_EL = real(P_etoile_A)-(imag(P_etoile_A)/tand(Phi_Z_EL));
P_EL = real(P_etoile_A)-(imag(P_etoile_A)/tand(Phi_P_EL));

%Cree une sous fonction de transfert pour trouver le Ka
TF_Ka_EL = tf([1 -Z_EL], [1 -P_EL]) * TF_EL;

%Calcul du K_AvPh_AZ
K_AvPh_EL = 1/abs(evalfr(TF_Ka_EL, P_etoile_A));

%Nouvelle fonction de transfert d'avance de phase
TF_AvPh_EL = TF_Ka_EL * K_AvPh_EL;

%% Calcul pour PI Elevation Télescope A
Diviser = 10;

        %Valider c'est de combien qu'on ajouter
        % figure
        % hold on
        % pzmap(TF_EL)
        % scatter(real(P_etoile_A), imag(P_etoile_A), 100, "square", 'black')

%Trouver les poles et zeros
Z_EL = real(P_etoile_A) / Diviser;

%Cree une sous fonction de transfert pour trouver le Kpi
TF_Kpi_EL = tf([1 -Z_EL], [1 0]) * TF_AvPh_EL;

%Calcul du K_PI_EL
K_PI_EL = 1/abs(evalfr(TF_Kpi_EL, P_etoile_A));

%Nouvelle fonction de transfert de PI
TF_Finale_EL = TF_Kpi_EL * K_PI_EL;
TF_Finale_BF_EL = feedback(TF_Finale_EL, 1);

%% Demande pour rapport
%Rlocus du système avec les P desirer
% figure
% hold on
% rlocus(TF_EL, "red")
% scatter(real(P_etoile_A), imag(P_etoile_A), 50, "*", 'black')
% scatter(real(P_etoile_A), imag(-P_etoile_A), 50, "*", 'black')
% xlim([-1400, 300])
% ylim([-200, 200])

%Réponse à l'échelon unitaire
% figure;
% step(feedback(TF_Finale_EL,1));

%Erreur à une parabole uniaire du système
temps = [0:0.1:2];
Parabole = [0:0.1:2].^2;
y = lsim(TF_Finale_BF_EL, Parabole, temps);

% figure
% hold on
% box on
% plot(temps, y', "blue")
% plot(temps, Parabole-y', "red")
% plot(temps, Parabole, "black")
% legend(["Réponse", "Erreur", "Parabole"]);
% title("Erreur à une rampe unitaire");
% ylabel("Amplitude");
% xlabel("Time (secondes)")

%Diagramme de bode du système
% figure
% margin(TF_Finale_EL)

%Erreur sur la trajectoire
% figure
% hold on
% box on
% y = lsim(TF_Finale_BF_EL, utrk, ttrk);
% %plot(ttrk, y, "blue")
% plot(ttrk, utrk-y, "red")
% %plot(ttrk, utrk, "black")
% %legend(["Réponse", "Erreur", "Trajectoire"]);
% legend(["Erreur"])
% title("Erreur sur la trajectoire");
% ylabel("Amplitude");
% xlabel("Time (secondes)")

%% Validation système AZ
%On vérifie Mp < 30%    Tr<0.25sec      ts< 1.20sec
stepinfo(TF_Finale_BF_EL, RiseTimeThreshold=[0 1]);

%on vérifie GM > 10 dB      RM > 0.09s
[Gm, Pm, wcg, wcp] = margin(TF_Finale_EL);
Pm;
Gm = 20*log10(Gm);
Rm = (Pm/wcp)*(pi/180);

Z_EL;
P_EL;
 
K_AvPh_EL;
TF_AvPh_EL;

% Kvel_EL;
% K_etoile_EL;

% P_RePh_EL;
% Z_RePh_EL;

% K_RePh_EL;
% TF_RePh_AvPh_EL;

%disp("Temps Erreur Parabole");
lsiminfo(Parabole-y', temps);

                    %Effacer les non utiliser
                    clear Gm Pm wcg wcp
                    %Clear les variables
                    clear Kvel_EL Kvel_etoile_EL K_etoile_EL Z_RePh_EL P_RePh_EL TF_Ka_EL
                    clear P_etoile_A Ajout_Omega_n Ajout_Omega_a
                    %Clear les variables
                    clear TF_Ka_EL P_AZ Z_EL Phi_P_EL Phi_Z_EL Angle_EL Delta_Phi_EL Alpha_EL


%Assurer la fin du document
disp("Hello World")