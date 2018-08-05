clear all; clc;

Pt = 1; n = 4; B = 10e6; Rb = 10e6; N0_db = -204; fc = 3e9; m = 2;
d = 10:1e3; c = 3e8; d0= 1;

lambda = c/fc;
Pt_db = 10*log10(Pt);
N0 = 10^(N0_db/10);

for i=1:length(d)
    Pl_el = (4*pi*d0/lambda)^2;
    Pl_el_db = 10*log10(Pl_el);
    Pl_ld_db = Pl_el_db + 10*n*log10(d(i)/d0);
    
    Pr_db = Pt_db - Pl_ld_db;
    Pr = 10^(Pr_db/10);
    
    SNR_avg(i) = Pr/(N0*B);
    
    C(i) = B*log10(1+SNR_avg(i));
    C(i) = Rb;
    SNR_min = 2^(C(i)/B)-1;
    
    P_rayl_aprox(i) = SNR_min/SNR_avg(i);
    P_rayl_exato(i) = 1-exp(-SNR_min/SNR_avg(i));
    P_naka_exato(i) = gammainc(SNR_min*m/SNR_avg(i),m,'lower')/gamma(m);
    P_naka_aprox(i) = (SNR_min*m/SNR_avg(i))^m/gamma(m+1);
end
figure;
semilogy(d, P_naka_exato,'k-','LineWidth',2);
hold on;
semilogy(d, P_naka_aprox,'k--','LineWidth',2);
semilogy(d, P_rayl_exato,'b-','LineWidth',2);
semilogy(d, P_rayl_aprox,'b--','LineWidth',2);

grid on;
title('P_{out} X d (Rayleigh e Nakagami)');
legend('Nakagami m = 2 exato','Nakagami m = 2 aprox','Rayleigh exato','Rayleigh aprox');
ylabel('P_{out}');
xlabel('Distância (m)');
%% Q3. Figura 4.2 Goldsmith

SNR_avg = 10^(20/10);
P_rayl = 1e-4:1e-3:1;
SNR_min = -log(1-P_rayl)*SNR_avg;
CB = log2(1+SNR_min);

figure;
semilogx(P_rayl,CB,'LineWidth',2);
grid on;
legend('C/B');
title('Capacidade Normalizada X P_{out}');
ylabel('C/B'); xlabel('P_{out}');

% Figura 4.2 Goldsmith c/ C0
figure;
SNR_min = -log(1-P_rayl).*SNR_avg;
CB = (1-P_rayl).*log2(1+SNR_min);

semilogx(P_rayl,CB,'LineWidth',2);
grid on;
%axis([1e-4 1 0 4]);
title('Capacidade Normalizada X P_{out}');
ylabel('C/B'); xlabel('P_{out}');
legend('C_{0}/B');

%% Q4. Capacidade de outage X SNR_avg
clear all; close all; clc;

SNR_avg_db = 0:300;
SNR_avg = 10.^(SNR_avg_db/10);
P_rayl = [1e-5 1e-4 1e-3 1e-2 1e-1];

legendCell = cell(length(P_rayl), 1);

for i=1:length(P_rayl)
    SNR_min = -log(1-P_rayl(i)).*SNR_avg;
    CB = log2(1+SNR_min);
    
    semilogy(SNR_avg_db,CB,'LineWidth',2);
    hold on;
    legendCell{i,1} = ['P_{out} = ' num2str(P_rayl(i))];
end

grid on;
legend(legendCell{:,1});
title('Capacidade Normalizada (C/B) X SNR média');
ylabel('C/B');
xlabel('SNR média');