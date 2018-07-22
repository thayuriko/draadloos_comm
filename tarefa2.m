%% Questão 1 - 8
close all; clear all; clc;

c = 3e8; fc = 2e9; d = 1e3; N0 = 1e-20; d0 = 1; n = 3;
fs = 44.1e3; l = 16; Rb1 = 1e6; Rb2 = fs*l;

M = [2, 4, 8, 16, 32];
B = 4.264e6;
%Pr/N0 = EbN0*Rb;
lambda = c/fc;
Rb = 10*(Rb1+Rb2);
BER = 1e-5;
Rs = B*2;

%BER = 2/log2(M(i))*qfunc(sqrt(2*EbN0*log2(M(i)))*sin(pi/M(i)));
%BER = 2*(M(i)-1)/(M(i)*log2(M(i)))*qfunc(sqrt(3*EbN0*log2(M(i))/(M(i)-1)));

for i=1:length(M)
    Rs(i) = Rb/log2(M(i));
    
    if Rs(i) == B
        % M-PSK
        EbN0_psk = ((qfuncinv(BER*log2(M(i))/2)/sin(pi/M(i)))^2)/(2*log2(M(i)));
        EbN0_db_psk(i) = 10*log10(EbN0_psk);
        Pr = EbN0_psk*Rb*N0;
        
        Pt = Pr*(16*pi^2*d^2)/lambda^2; %espaco livre
        Pr_db = 10*log10(Pr);
        Pt_dbm(i,1) = 10*log10(Pt) + 30;

        Pl_el = (16*pi^2*d0^2)/lambda^2;
        Pl_el_db = 10*log10(Pl_el);
        Pl_ld_db = Pl_el_db + 10*n*log10(d/d0);
        Pt_db = Pr_db + Pl_ld_db; %log distância
        Pt_dbm(i,2) = Pt_db + 30;
        
        
        % M-QAM
        EbN0_qam = BER*(M(i)*log2(M(i)))/(2*(M(i)-1));
        EbN0_qam = (M(i)-1)*qfuncinv(EbN0_qam)^2/(3*log2(M(i)));
        EbN0_db_qam(i) = 10*log10(EbN0_qam);
        
        Pr = EbN0_qam*Rb*N0;
        
        Pt = Pr*(16*pi^2*d^2)/lambda^2; %espaco livre
        Pr_db = 10*log10(Pr);
        Pt_dbm(i,3) = 10*log10(Pt) + 30;

        Pl_el = (16*pi^2*d0^2)/lambda^2;
        Pl_el_db = 10*log10(Pl_el);
        Pl_ld_db = Pl_el_db + 10*n*log10(d/d0);
        Pt_db = Pr_db + Pl_ld_db; %log distância
        Pt_dbm(i,4) = Pt_db + 30;
        
        fprintf('Q1. M = %s\n', num2str(M(i)));
        fprintf('Q2. Pt_dbm = %s dBm\n', num2str(Pt_dbm(i,1)));
        fprintf('Q3. M = %s\n', num2str(M(i)));
        fprintf('Q4. Pt_dbm = %s dBm\n', num2str(Pt_dbm(i,3)));
        fprintf('Q5. M = %s\n', num2str(M(i)));
        fprintf('Q6. Pt_dbm = %s dBm\n', num2str(Pt_dbm(i,2)));
        fprintf('Q7. M = %s\n', num2str(M(i)));
        fprintf('Q8. Pt_dbm = %s dBm\n', num2str(Pt_dbm(i,4)));
    else
        Pt_dbm(i,1) = 0; %psk espaço livre
        Pt_dbm(i,2) = 0; %psk log distância
        Pt_dbm(i,3) = 0; %qam espaço livre
        Pt_dbm(i,4) = 0; %qam log distância
    end 
end

%% Questão 9 - 14
clear all;

%Pl_db = -50 + 10*log10(fc) + 30*log10(d);
fc = 1e9; N0_dbm = -174; B = 1e6; EsN0_db = 10;
c = 3e8;

N0_db = N0_dbm - 30;
N0 = 10^(N0_db/10);
EsN0 = 10^(EsN0_db/10);

Rs = B;

%Q.9
Pr = EsN0*Rs*N0;
Pr_db = 10*log10(Pr);
fprintf('Q9. Pr_min = %s dB\n', num2str(Pr_db));

%Q.10
Pt = 1;
Pt_db = 10*log10(Pt);

%Pl_db = Pt_db - Pr_db = -50 + 10*log10(fc) + 30*log10(d);
d = 10^(((Pt_db - Pr_db) + 50 - 10*log10(fc))/30);
fprintf('Q10. d = %s m\n', num2str(d));

%Q.11
d = 2e3;
Pl_db = -50 + 10*log10(fc) + 30*log10(d);
Pt_db = Pl_db + Pr_db;
fprintf('Q11. Pt = %s dB\n', num2str(Pt_db));

%Q.12
Rs = 100e3;
Pr = EsN0*Rs*N0;
Pr_db = 10*log10(Pr);
Pt = 1;
Pt_db = 10*log10(Pt);
d_new = 10^(((Pt_db - Pr_db) + 50 - 10*log10(fc))/30);
fprintf('Q12. d = %s vezes\n', num2str(d_new/d));

%Q.13
fprintf('Q13. 0.5\n');

%Q.14
sigma = 4;
%qfunc((gama - mi)/sigma) = 1e-3;
Pt = qfuncinv(1e-3)*sigma;
fprintf('Q14. Pt = %s dB\n',num2str(Pt));

%% Questão 15 - 17
clear all;

%Pr_db=Pt_db+(Gt_db+Gr_db)-Pl_db; -> Pt_db=Pr_db-(Gt_db+Gr_db)+Pl_db;

% Cálculo do sistema A:
Gt_db = 0; Pl_db = 110; Rs = 10e3; Gr_db = 2.3; N0_db = -195;
margem_db = 10; BER = 1e-3;

EbN0 = ((qfuncinv(BER*log2(8)/2)/sin(pi/8))^2)/(2*log2(8));

EbN0_db = 10*log10(EbN0);
EbN0mg_db = EbN0_db + margem_db;
EbN0mg = 10^(EbN0mg_db/10);
N0 = 10^(N0_db/10);
Pr = EbN0mg*Rs*N0;
Pr_db = 10*log10(Pr);
PtA_db=Pr_db-(Gt_db+Gr_db)+Pl_db;
fprintf('Q15. PtA = %s dB\n', num2str(PtA_db));

% Cálculo do sistema B:
Gt_db = 0; Pl_db = 110; Rs = 10e5; Gr_db = 2.3; N0_db = -195;
margem_db = 10; BER = 1e-2;

EbN0 = ((qfuncinv(BER*log2(4)/2)/sin(pi/4))^2)/(2*log2(4));

EbN0_db = 10*log10(EbN0);
EbN0mg_db = EbN0_db + margem_db;

EbN0mg = 10^(EbN0mg_db/10);
N0 = 10^(N0_db/10);
Pr = EbN0mg*Rs*N0;
Pr_db = 10*log10(Pr);
PtB_db=Pr_db-(Gt_db+Gr_db)+Pl_db;
fprintf('Q16. PtB = %s dB\n', num2str(PtB_db));

%% Questão 18
clear all;

Pt = 4e-3; Gt_db = 5; Gr_db = 8; Pl_db = 100; margem_db = 5; 
N0_db = -204; BER = 0.01e-2; B = 5e6; alpha = 0.25; Bmax = 0;

Pt_db = 10*log10(Pt);
Pl = 10^(Pl_db/10);
N0 = 10^(N0_db/10);
Pr_db = Pt_db + (Gt_db + Gr_db) - Pl_db - margem_db;
Pr = 10^(Pr_db/10);

M = [2, 4, 8, 16, 32];

%BER = 2/log2(M(i))*qfunc(sqrt(2*EbN0*log2(M(i)))*sin(pi/M(i)));
for i=1:length(M)
    EbN0 = ((qfuncinv(BER*log2(M(i))/2)/sin(pi/M(i)))^2)/(2*log2(M(i)));
    
    Rb1 = Pr/(N0*EbN0);
    Rs = Rb1/log2(M(i));
    Bcl(i) = Rs*(1+alpha);
    
    Rs = B/(1+alpha);
    Rb2 = Rs*log2(M(i));
    
    if B >= Bcl(i)
        if Bcl(i) > Bmax
            Bmax_id = i;
            Bmax = Bcl(i);
        end
        
        Rb_kbps(i,1) = Rb1/1e3;
        Rb_kbps(i,2) = Rb2/1e3;
    end
end

fprintf('Q18. Rb_max = %s kbps\n', num2str(max(Rb_kbps(Bmax_id,:))));

%% Questão 19
clear all;

Pt = 100e-3; Pr_dbm = -105; d0 = 50; d = 500; n = 4.2;
Gt_db = -7; Gr_db = 9; fc = 900e6; c = 3e8;

lambda = c/fc;
Pr_db = Pr_dbm - 30;
Pt_db = 10*log10(Pt);

%Pr_db = Pt_db + G_db - Pl_db;

Pl_el_db = 20*log10(4*pi*d0/lambda);
Pl_ld_db = Pl_el_db + 10*n*log10(d/d0);
Pr_db2 = Pt_db + (Gr_db + Gt_db) - Pl_ld_db;

margem_db = Pr_db2 - Pr_db;

fprintf('Q19. margem_db = %s dB\n', num2str(margem_db));