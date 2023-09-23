% Filtre de kalman pour l’estimation de position

clear variables;
close all;
clc;

load('mesurestrajKalm2.mat');
% mesures polaires
% sigmesurayon
% simesuang
% Tmesu

Z(:,1) = mesures(:,1) .* cos(mesures(:,2));
Z(:,2) = mesures(:,1) .* sin(mesures(:,2));

N = size(Z,1);

% echantillonage de 0.1s
T = 0.1;

g = 9.81;

sigma_a = 0.05;
mu = 0.01;
beta = exp(-mu);
sigma_w = sigma_a * sqrt(1-exp(-2*mu));

% FILTRE DE KALMAN
% initialization
X0 = [0;0;0;0;0;g*(beta - 1)];
Pi = zeros(6,6,N);
P0 = [10^6 0 0 0 0 0;
        0 10^6 0 0 0 0;
        0 0 10^6 0 0 0;
        0 0 0 10^6 0 0;
        0 0 0 0 10^6 0;
        0 0 0 0 0 10^6];
F = [1 0 T 0 T^2/2 0 ;
    0 1 0 T 0 T^2/2 ;
    0 0 1 0 T 0 ;
    0 0 0 1 0 T ;
    0 0 0 0 beta 0 ;
    0 0 0 0 0 beta ];
H = [1 0 0 0 0 0 ;
    0 1 0 0 0 0 ];
W = T * [0 0 0 0 0 0 ;
        0 0 0 0 0 0 ;
        0 0 0 0 0 0 ;
        0 0 0 0 0 0 ;
        0 0 0 0 sigma_w^2 0 ;
        0 0 0 0 0 sigma_w^2];

X_est = zeros(6,N);
var = zeros(1,N);
erreur = zeros(1,N);

% ---------
% main loop
% ---------
for i = 1:N-1
    % prediction
    if i > 1
        X_est(:,i) = F*X_est(:,i-1);
        Pi(:,:,i) = F*Pi(:,:,i-1)*F' + W;
    else
        X_est(:,1) = F*X0;
        Pi(:,:,1) = F*P0*F' + W;
    end
    erreur(1,i) = sum(sqrt((H*X_est(:,i)-Z(i,:)').*(H*X_est(:,i)-Z(i,:)')));
    % correction
    % gain de Kalman
    D = mesures(i,1);
    alpha = mesures(i,2);
    J = [cos(alpha) -D*sin(alpha) ;
        sin(alpha) D*cos(alpha)];
    V = [sigmesurayon^2 0 ;
        0 simesuang^2];
    Vn = J*V*J';
    Ki = Pi(:,:,i)*H'*inv(H*Pi(:,:,i)*H'+Vn);
    % mise a jour
    X_est(:,i) = X_est(:,i) + Ki*(Z(i+1,:)'-H*X_est(:,i));
    Pi(:,:,i) = (eye(6)-Ki*H)*Pi(:,:,i);
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
plot(Tmesu(1:end-1),Z(1:end-1,1),Tmesu(1:end-1),X_est(1,1:end-1),'r')
legend('Donné réel','Donné estimé')

figure(4)
plot(Tmesu(1:end-1),Z(1:end-1,2),Tmesu(1:end-1),X_est(2,1:end-1),'r')
legend('Donné réel','Donné estimé')