clear all; close all; clc;

Pt_dbm = 12;

load('cc1350_1.mat');
load('cc1350_2.mat');

pathloss(1,:) = dados1_m(:,3);
pathloss(2,:) = dados10_m(:,3);
pathloss(3,:) = dados20_m(1:1000,3);
%pathloss(3,:) = repmat(mean(dados20_m(:,3)),1,1000);
pathloss(4,:) = dados30_m(:,3);
pathloss(5,:) = dados40_m(:,3);

fading_dbm(1,:) = dados(:,3);
fading_dbm(2,:) = dados20_m(:,3);

syms x f(x);
d = [1 10 20 30 40];
Pr_dbm = [];

clear dados1_m dados10_m dados20_m dados30_m dados40_m dados;
f(x) = 0;

for i=1:size(pathloss,1)
    Pr_dbm(i) = mean(pathloss(i,:));
    Ei(i) = Pr_dbm(1) - 10*x*log10(d(i));
    
    f(x) = f(x) + (Pr_dbm(i) - Ei(i))^2;
end

df = diff(f,x);
n = double(solve(df,x));

Pr_sim = Pr_dbm(1) - 10*n*log10(d);

figure;
plot(d,Pr_dbm, '-o', 'LineWidth', 2); hold on;
plot(d,Pr_sim, '-s', 'LineWidth', 2); hold off;
legend('Medido', 'Simulado');
title(['Perda de Percurso: Modelo Log-distância (n = ' num2str(n) ')']);
ylabel('RSSI (dB)'); xlabel('distância (m)');
grid on;

%{
1.a) Há uma divergência entre as curvas estimada e medida, no entanto, isso
ocorre pelo fato da curva estimada descrever o comportamento médio da perda
de percurso e não valores instantâneos de RSSI.

1.b) As medições foram realizadas em um ambiente do tipo espaço-livre,
condizente com o valor de n estimado, não podendo, portanto, ser utilizado
em outros ambientes. O coeficiente indica a velocidade com a qual a
perda de percurso aumenta em relação à distância entre Rx e Tx, dependendo
predominantemente do ambiente de propagação.
%}

%%
figure;
plot(fading_dbm(2,:)); hold on;
plot(fading_dbm(1,:)); hold off;
legend('Medida 1', 'Medida 2');
title('RSSI instantânea medida para 10.000 amostras');
grid on; xlabel('Amostras'); ylabel('RSSI (dB)');

fading = 10.^(fading_dbm/10);
h1 = sqrt(fading(1,:)/mean(fading(1,:)));
h2 = sqrt(fading(2,:)/mean(fading(2,:)));
dfittool;

ylabel('pdf(h)'); xlabel('h');
title('PDF da Distribuição Nakagami-m'); grid on;
legend('Medidas', ['Nakagami m = ' num2str(pd_1.mu)]);