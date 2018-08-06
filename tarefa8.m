%% Q1. P_out em função da distância utilizando MRC
% github.com/thayuriko/draadloos_comm
clear all; close all; clc;
addpath('img')

Pt = 1; n = 4; B = 10e6; Rb = 10e6; N0_db = -204; fc = 3e9; m = 2;
d = 1:5e2; c = 3e8; d0 = 1;
nRx = [1 2 4 8];

lambda = c/fc;
Pt_db = 10*log10(Pt);
N0 = 10^(N0_db/10);

Pl_el = (4*pi*d0/lambda)^2;
Pl_el_db = 10*log10(Pl_el);

for j=1:length(nRx)
    for i=1:length(d)
        Pl_ld_db = Pl_el_db + 10*n*log10(d(i)/d0);

        Pr_db = Pt_db - Pl_ld_db;
        Pr = 10^(Pr_db/10);

        SNR_avg(i) = Pr/(N0*B);
        SNR_min = 2^(Rb/B)-1;

        P_mrc(j,i) = gammainc(SNR_min*nRx(j)/SNR_avg(i),nRx(j),'lower');
    end
end

plot_mrc = figure('units','normalized','outerposition',[0 0 0.5 1]);
semilogy(d, P_mrc(1,:),'-','LineWidth',2); hold on;
semilogy(d, P_mrc(2,:),'-','LineWidth',2); 
semilogy(d, P_mrc(3,:),'-','LineWidth',2); 
semilogy(d, P_mrc(4,:),'-','LineWidth',2); hold off;
grid on; axis([1 d(end) 1e-8 1]);
title('Probabilidade de Outage no canal Rayleigh com MRC');
legend('nRx = 1','nRx = 2','nRx = 4','nRx = 8');
ylabel('P_{out}'); xlabel('Distância (m)');

print(plot_mrc,'img/t8_q1','-dpng');
%% Q2. P_out em função da SNR utilizando MRC e SC
% github.com/thayuriko/draadloos_comm
clear all; close all; clc;
addpath('img')

B = 10e6; Rb = 10e6;
nRx = [1 2 4 8];
SNR_avg_db = -10:60;
SNR_avg = 10.^(SNR_avg_db./10);

for j=1:length(nRx)
    for i=1:length(SNR_avg_db)
        SNR_min = 2^(Rb/B)-1;

        P_mrc(j,i) = gammainc(SNR_min*nRx(j)/SNR_avg(i),nRx(j),'lower');
        P_sc(j,i) = (1-exp(-SNR_min/SNR_avg(i)))^nRx(j);
    end
end

plot = figure('units','normalized','outerposition',[0 0 0.5 1]);
semilogy(SNR_avg_db, P_mrc(1,:),'-','LineWidth',2); hold on;
semilogy(SNR_avg_db, P_mrc(2,:),'-','LineWidth',2); 
semilogy(SNR_avg_db, P_mrc(3,:),'-','LineWidth',2); 
semilogy(SNR_avg_db, P_mrc(4,:),'-','LineWidth',2); 
semilogy(SNR_avg_db, P_sc(1,:),'--','LineWidth',2);
semilogy(SNR_avg_db, P_sc(2,:),'--','LineWidth',2); 
semilogy(SNR_avg_db, P_sc(3,:),'--','LineWidth',2); 
semilogy(SNR_avg_db, P_sc(4,:),'--','LineWidth',2); hold off;
grid on; axis([SNR_avg_db(1) SNR_avg_db(end) 1e-15 1]);
title('Probabilidade de Outage no canal Rayleigh com MRC e SC');
legend('MRC nRx = 1','MRC nRx = 2','MRC nRx = 4','MRC nRx = 8', ...
    'SC nRx = 1','SC nRx = 2','SC nRx = 4','SC nRx = 8');
ylabel('P_{out}'); xlabel('SNR (dB)');

print(plot,'img/t8_q2','-dpng');