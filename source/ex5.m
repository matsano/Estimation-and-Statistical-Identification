% Filtre de kalman pour l’estimation de position et vitesse

clear variables;
close all;
clc;

load('mesurestrajKalm3.mat');
% mesures polaires
% sigmesurayon
% simesuang
% Tmesu

N = size(mesures,1);

% echantillonage de 0.1s
T = 0.1;

g = 9.81;

sigma_a = 0.05;
mu = 0.01;
beta = exp(-mu);
sigma_w = sigma_a * sqrt(1-exp(-2*mu));

% FILTRE DE KALMAN
% initialisation
x_0 = mesures(1,1) * cos(mesures(1,6));
y_0 = mesures(1,1) * sin(mesures(1,6));
X0 = [x_0/2;y_0/2;0;0;0;-g;0;0];
%X0 = [0;0;0;0;0;-g;0;0];
Pi = zeros(8,8,N);
P0 = 10^6 * eye(8);
F = [1 0 T 0 T^2/2 0 0 0 ;
    0 1 0 T 0 T^2/2 0 0 ;
    0 0 1 0 T 0 0 0 ;
    0 0 0 1 0 T 0 0 ;
    0 0 0 0 beta 0 0 0 ;
    0 0 0 0 0 beta 0 0 ;
    0 0 0 0 0 0 0 0 ;
    0 0 0 0 0 0 0 0 ];
H = [0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 1];
% TODO
W = T * [0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 sigma_w^2 0 0 0;
        0 0 0 0 0 sigma_w^2 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0];
V = [sigmesurayon^2 0 ;
    0 simesuang^2];
X_est = zeros(8,N);
var = zeros(1,N);
erreur = zeros(1,N);

X_est_polar = zeros(200,2);
eqm = zeros(5,1);
goodMesure = zeros(200,2);

% ---------
% main loop
% ---------
for i = 1:N
    % prediction
    if i > 1
        X_est(:,i) = F*X_est(:,i-1);
        Pi(:,:,i) = F*Pi(:,:,i-1)*F' + W;
    else
        X_est(:,1) = F*X0;
        Pi(:,:,1) = F*P0*F' + W;
    end

    x = X_est(1,i);
    y = X_est(2,i);

    % linearization matrice H

    Hk = [x/sqrt(x^2+y^2) y/sqrt(x^2+y^2)
        -y/(x^2+y^2) x/(x^2+y^2)];
    H(1:2,1:2) = Hk;

    X_est_polar(i,:) = [sqrt(x^2+y^2); atan(y/x)];
    
    % mise a jour mk_x et mk_y
    X_est(7:8,i) = -H*X_est(:,i) + [sqrt(x^2+y^2); atan(y/x)];
    
    % NNSF - selection des mesures les plus probables
    xk = X_est_polar(i,1);
    yk = X_est_polar(i,2);
    for j = 1:5
    eqm(j) = sum(sqrt(([xk;yk]-[mesures(i,j);mesures(i,j+5)]).^2));
    end
    [erreur(1,i), j_best] = min(eqm);
    goodMesure(i,:) = [mesures(i,j_best) mesures(i,j_best+5)];
    
    % correction
    % gain de Kalman
    Ki = Pi(:,:,i)*H'*inv(H*Pi(:,:,i)*H'+V);
    % mise a jour
    X_est(:,i) = X_est(:,i) + Ki*(goodMesure(i,:)'-H*X_est(:,i));
    Pi(:,:,i) = (eye(8)-Ki*H)*Pi(:,:,i);
    var(1,i) = trace(Pi(:,:,i));
end

% plot
n = 1:N;

figure(1)
plot(n,var);
moyenne = mean(var(end-10:end));
str = strcat('AVG variance: ',num2str(moyenne));
text(130,500000,str);

figure(2)
plot(n,erreur)
moyenne = mean(erreur(end-10:end));
str = strcat('AVG error: ',num2str(moyenne));
text(130,1000,str);

figure(3)
plot(Tmesu(1:end-1),goodMesure(1:end-1,1),Tmesu(1:end-1),X_est_polar(1:end-1,1),'r')
legend('Donné réel','Donné estimé')

figure(4)
plot(Tmesu(1:end-1),goodMesure(1:end-1,2),Tmesu(1:end-1),X_est_polar(1:end-1,2),'r')
legend('Donné réel','Donné estimé')
