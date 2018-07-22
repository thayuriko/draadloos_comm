%% Questão 1
close all; clear all; clc;

c = 3e8; fc = 900e6; d = 10;
Gt_db = 1; Gr_db = 1; Pr = 10e-6;

Gt = 10^(Gt_db/10);
Gr = 10^(Gr_db/10);
lambda = c/fc;

Pl = 16*pi^2*d^2/lambda^2;
Pl_db = 10*log10(Pl);
Pr_db = 10*log10(Pr);

Pt_db = Pr_db + Pl_db;

fprintf('Questão 1.\nA potência de transmissão é de %sdB\n\n', num2str(Pt_db));

%% Questão 2
clear all;

c = 3e8; fc = 5e9; d = 10;
Gt_db = 1; Gr_db = 1; Pr = 10e-6;

Gt = 10^(Gt_db/10);
Gr = 10^(Gr_db/10);
lambda = c/fc;

Pl = 16*pi^2*d^2/lambda^2;
Pl_db = 10*log10(Pl);
Pr_db = 10*log10(Pr);

Pt_db = Pr_db + Pl_db;

fprintf('Questão 2.\nA potência de transmissão é de %sdB\n\n', num2str(Pt_db));

%% Questão 3
clear all;

c = 3e8; fc = 1800e6; d0 = 1e3; d = 1e3:1e3:20e3;
Pr0 = 1e-6; Gt_db = 1; Gr_db = 1; n = [3 4];
ht = 40; hr = 3;

lambda = c/fc;
Gt = 10^(Gt_db/10);
Gr = 10^(Gr_db/10);
Pr0_db = 10*log10(Pr0);
ahr_db = 3.2*(log10(11.75*hr))^2-4.97;
Cm = 3;

for i=1:length(d)
    Pr(i,1) = Pr0*(d0/d(i))^2; % Espaço livre
    Pt = Pr0*(16*(pi^2)*(d0^2))/(Gt*Gr*lambda);
    Pt_db = 10*log10(Pt);
    
    Pr(i,2) = Pr0*(d0/d(i))^n(1); % Log distância
    Pr(i,3) = Pr0*(d0/d(i))^n(2); % Log distância
    
    Pl0_db = 46.3+33.9*log10(fc/1e6)-13.82*log10(ht)-ahr_db+(44.9-6.55*log10(ht))*log10(d0/1e3)+Cm; %Hata
    Pt_db = Pr0_db + Pl0_db;
    
    Pl_db = 46.3+33.9*log10(fc/1e6)-13.82*log10(ht)-ahr_db+(44.9-6.55*log10(ht))*log10(d(i)/1e3)+Cm; %Hata
    Pl = 10^(Pl_db/10);
    Pr_db = Pt_db - Pl_db;
    
    Pr(i,4) = 10^(Pr_db/10); %Okumura-Hata estendido
end

Pr_db = 10*log10(Pr);

figure;
loglog(d,Pr(:,1),'*-'); hold on;
loglog(d,Pr(:,2),'d-');
loglog(d,Pr(:,3),'s-');
loglog(d,Pr(:,4),'x-');
axis([1e3 20e3 1e-12 1e-6]);
grid on;

legend('Espaço Livre', 'Log Distância n = 3', 'Log Distância n = 4', 'Okumura-Hata estendido');
title('Desempenho de diferentes modelos de perda de percurso');
xlabel('Distância (km)');
ylabel('Potência recebida (W)');

%% Questão 4
clear all;

c = 3e8;
d0 = 1; Pt = 1e-3; fc = 900e6; d = 200;
table = [10 -70; 20 -75; 50 -90; 100 -110; 300 -125];

lambda = c/fc;

Pl0_db = -20*log10(lambda/(4*pi));
%(a+bx)^2 = a^2+2abx+x^2;
a = 0; b = 0; ab = 0;

for i=1:length(table)
    a = (table(i,2)-Pl0_db)^2 + a;
    b = (10*log10(table(i,1)))^2 + b;
    ab = 2*(table(i,2)+Pl0_db)*(10*log10(table(i,1))) + ab;
end

n = -ab/(2*b);

Pr0 = Pt*lambda^2/(16*pi^2*d0^2);
Pr0_db = 10*log10(Pr0);

Pr_db = Pr0_db - 10*n*log10(d/d0);
Pr_dbm = Pr_db + 30;

fprintf('Questão 4.\nA potência recebida a 150m é de %sdBm\n\n', num2str(Pr_dbm));

%% Questão 5

for i=10:10:300
    Pr_db(i/10) = Pr0_db - 10*n*log10(i/d0);
    Pr(i/10) = 10^(Pr_db(i/10)/10);
end

figure;
loglog([10:10:300],Pr,'*-');
grid on;
axis([10 300 1e-16 1e-9]);
title('Modelo de perda de percurso log-normal');
xlabel('Distância (m)');
ylabel('Potência recebida (W)');