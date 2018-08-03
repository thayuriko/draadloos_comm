%% Q1
clear all; close all; clc;

Pt = 1; n = 4; B = 10e6; Rb = 10e6; N0_db = -204; fc = 3e9; m = 2;
d = 10:1e3; c = 3e8; d0= 1;
nRx = [1 2 4 8];

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
    
    for j=1:length(nRx)
        P_rayl(j,i) = gammainc(nRx(j)*SNR_min/SNR_avg(i),nRx(j),'lower')/gamma(nRx(j));
    end
end

figure;
semilogy(d, P_rayl(1,:),'-','LineWidth',2); hold on;
semilogy(d, P_rayl(2,:),'-','LineWidth',2); 
semilogy(d, P_rayl(3,:),'-','LineWidth',2); 
semilogy(d, P_rayl(4,:),'-','LineWidth',2); hold off;

%{
SNR_avg_db = 10*log10(SNR_avg);
figure
semilogy(SNR_avg_db, P_rayl(1,:),'-','LineWidth',2); hold on;
semilogy(SNR_avg_db, P_rayl(2,:),'-','LineWidth',2);
semilogy(SNR_avg_db, P_rayl(3,:),'-','LineWidth',2);
grid on;
%}
grid on;
%axis([0 60 1e-5 1]);
title('P_{out} X d (Rayleigh + MRC)');
legend('Rayleigh 1 nRx','Rayleigh nRx 2','Rayleigh nRx 4');
ylabel('P_{out}');
xlabel('Distância (m)');
%% Q2
clear all; close all; clc;

rand('state',0); randn('state',0); bar = waitbar(0,'Setting up...');

EbN0_db = 0:40;
nBits = 1e7;
EbN0 = 10.^(EbN0_db/10);
[ber_awgn, ber_rayl] = deal(zeros(1,length(EbN0)));
nRx = [1 2 4 8];  %vetor com a qtde de antenas

m = rand(1,nBits)>0.5;
x = 2*m-1;      %BPSK

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait');
    
    for j=1:length(nRx)
        N0 = 1/EbN0(i);
        n = sqrt(0.5*N0)*(randn(nRx(j),nBits)+1i*randn(nRx(j),nBits));
        h = sqrt(0.5)*(randn(nRx(j),nBits)+1i*randn(nRx(j),nBits));

        sd = kron(ones(nRx(j),1),x);
        y = sd.*h + n;

        yHat =  sum(conj(h).*y,1)./sum(h.*conj(h),1); 
        err2 = sum((yHat>0)~= m);
        ber_rayl(j,i) = err2/nBits;
    end
    
    y0 = x + n(1,:);
    err0 = sum((y0>0)~= m);
    ber_awgn(i) = err0/nBits;
end

figure;
semilogy(EbN0_db, ber_rayl(1,:), '-o', 'Linewidth', 2); hold on;
semilogy(ber_rayl(2,:), '-d', 'Linewidth', 2);
semilogy(ber_rayl(3,:), '-s', 'Linewidth', 2);
semilogy(ber_awgn, '-*', 'Linewidth', 2);
grid on; axis([1 length(EbN0)-1 1e2/nBits 1]);
ylabel('BER'); xlabel('EbN0 (dB)');
legend('Rayleigh + MRC (nRx = 1)', 'Rayleigh + MRC (nRx = 2)', 'Rayleigh + MRC (nRx = 4)', 'AWGN');
title('BER da BPSK em um canal Rayleigh utilizando a técnica MRC');
close(bar);