clear all; close all; clc;

Pt_dbm = 12; Rb = 50e3;

load('cc1350_1.mat');
load('cc1350_2.mat');

pathloss(1,:) = dados1_m(:,3);
pathloss(2,:) = dados10_m(:,3);
pathloss(3,:) = dados20_m(1:1000,3);
%pathloss(3,:) = repmat(mean(dados20_m(:,3)),1,1000);
pathloss(4,:) = dados30_m(:,3);
pathloss(5,:) = dados40_m(:,3);

fading(1,:) = dados20_m(:,3);
fading(2,:) = dados20_m(:,3);

syms n f(x) x j(n);
d = [1 10 20 30 40];
Pr_dbm = [];

clear dados1_m dados10_m dados20_m dados30_m dados40_m dados;

for i=1:size(pathloss,1)
    Pr_dbm(i) = mean(pathloss(i,:));
    Ei(i) = Pr_dbm(1) - 10*n*log10(d(i)/d(1));
end

j(n) = sum(Pr_dbm - Ei)^2;
df = diff(j,n);
coef = double(solve(df,n));

Pr_sim = Pr_dbm(1) - 10*2.586*log10(d./d(1));

figure;
plot(d,Pr_dbm, '-o', 'LineWidth', 2); hold on;
plot(d,Pr_sim, '-s', 'LineWidth', 2); hold off;
legend('Medido', 'Simulado');
title(['Perda de Percurso: Modelo Log-distância (n = ' num2str(coef) ')']);
ylabel('RSSI (dB)'); xlabel('distância (m)');
grid on;