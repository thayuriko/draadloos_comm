%% Q2 - BER e FER do BPSK com códigos com desvanecimento Rayleigh (canais rápido e lento)

%---------------FAST FADING---------------%
clear all; close all; clc;

rand('state',0); randn('state',0);
bar = waitbar(0,'Setting up...');

EbN0_db = 1:30;
nBits = 10^6;
EbN0 = 10.^(EbN0_db/10);
ber_awgn = zeros(1,length(EbN0));
n = 15; k = 11;

disp('Desempenho Código de Hamming (15,11) p/ BPSK');

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait (fast fading)');
    N0 = 1/EbN0(i);
    m = rand(1,nBits)>0.5;  %geração de bits aleatórios
    [nBits_enc, err_enc] = deal(0);
    bpsk = 2*m-1;           %modulação bpsk sem codificação
    
    noise = sqrt(0.5*N0)*(randn(1,nBits)+1i*randn(1,nBits)); %awgn
    h = sqrt(0.5)*(randn(1,nBits)+1i*randn(1,nBits)); %rayleigh
    
    y = bpsk.*h + noise; %sinal transmitido
    
    err = sum((y./h>0)~= m);
    
    for j=1:k:(length(m)-rem(length(m),k))
        frame = m(j:j+k-1); %formação dos blocos
        frame_enc = encode(frame,n,k,'hamming/binary'); %codificação
        bpsk_enc = 2*frame_enc-1;   %bloco codificado e modulado
        
        noise = sqrt(0.5*N0)*(randn(1,n)+1i*randn(1,n));
        h = sqrt(0.5)*(randn(1,n)+1i*randn(1,n));
        
        y = bpsk_enc.*h + noise; %sinal transmitido
        
        y_dem = y./h>0; %demodulado
        y_dec = decode(y_dem,n,k,'hamming/binary'); %decodificado
        err_enc = sum(y_dec ~= frame) + err_enc;    %erros acumulados
        nBits_enc = nBits_enc + 1;
    end
    
    ber(i) = err/nBits;
    ber_enc(i) = err_enc/(nBits_enc*k);
end

fer = 1-(1-ber).^n;
fer_enc = 1-(1-ber_enc).^n;

plot_fast = figure('NumberTitle', 'off','units','normalized','outerposition',[0 0 1 1],'Name', 'Fast Fading');
subplot(1,2,1);
semilogy(EbN0_db,ber_enc,'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,ber,'-d','Linewidth', 2);
%title('BER x EbN0');
legend('Codificado - Hamming (15,11)', 'Não-codificado','location', 'best');
ylabel('BER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

subplot(1,2,2);
semilogy(EbN0_db,fer_enc,'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,fer,'-d','Linewidth', 2);
legend('Codificado - Hamming (15,11)', 'Não-codificado','location', 'best');
%title('FER x EbN0');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);
close(bar);

a = axes;
t1 = title('Performance BPSK com desvanecimento Rayleigh em Canal Rápido','FontSize',16);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');


%---------------SLOW FADING---------------%
bar = waitbar(0,'Setting up...');
ber_awgn = zeros(1,length(EbN0));

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait (slow fading)');
    N0 =  1/EbN0(i);
    m = rand(1,nBits)>0.5;  %geração de bits aleatórios
    [nBits_enc, err_enc, err] = deal(0);
    
    for j=1:k:(length(m)-rem(length(m),k))
        frame = m(j:j+k-1);                                 %formação dos blocos
        
        % Rayleigh não-codificado
        bpsk = 2*frame-1;                                   %bpsk
        noise = sqrt(0.5*N0)*(randn(1,k)+1i*randn(1,k));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,k); %rayleigh canal semi-estático
        y = bpsk.*h + noise;                                %transmitido
        err = sum((y./h>0)~= frame) + err;                  %soma de bits errados
        
        % Codificação Hamming (15,11)
        frame_enc = encode(frame,n,k,'hamming/binary');     %codificação
        bpsk_enc = 2*frame_enc-1;                           %bpsk
        noise = sqrt(0.5*N0)*(randn(1,n)+1i*randn(1,n));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,n); %rayleigh canal semi-estático
        y = bpsk_enc.*h + noise;                            %sinal transmitido
        y_dem = y./h>0;                                     %demodulado
        y_dec = decode(y_dem,n,k,'hamming/binary'); %decodificado
        err_enc = sum(y_dec ~= frame) + err_enc;    %erros acumulados
        nBits_enc = nBits_enc + 1;
    end
    
    ber(i) = err/nBits;
    ber_enc(i) = err_enc/(nBits_enc*k);
end

fer = 1-(1-ber).^n;
fer_enc = 1-(1-ber_enc).^n;

plot_slow = figure('NumberTitle', 'off','units','normalized','outerposition',[0 0 1 1],'Name', 'Slow Fading');
subplot(1,2,1);
semilogy(ber_enc,'-o','Linewidth', 2); hold on;
semilogy(ber,'-d','Linewidth', 2);
legend('Codificado - Hamming (15,11)', 'Não-codificado','location', 'best');
ylabel('BER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

subplot(1,2,2);
semilogy(EbN0_db,fer_enc,'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,fer,'-d','Linewidth', 2);
legend('Codificado - Hamming (15,11)', 'Não-codificado','location', 'best');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

a = axes;
t1 = title('Performance BPSK com desvanecimento Rayleigh em Canal Lento','FontSize',16);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');

print(plot_slow,'t6_q2_slow-Thais','-dpng')
print(plot_fast,'t6_q2_fast-Thais','-dpng')

close(bar);
%% Q4 - Simular a FER do BPSK com códigos considerando 0, 1 e 2 retransmissões (HARQ Tipo-I Simples)
%https://github.com/thayuriko/draadloos_ntw
clear all; close all; clc;

rand('state',0); randn('state',0);

EbN0_db = 0:30;
nBits = 10^6;
EbN0 = 10.^(EbN0_db/10);
ber_awgn = zeros(1,length(EbN0));
n = 15; k = 11;

bar = waitbar(0,'Setting up...');

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait');
    N0 =  1/EbN0(i);
    m = rand(1,nBits)>0.5;  %geração de bits aleatórios
    [err, nBits_enc] = deal(0);
    err_enc = zeros(1,3);
    
    for j=1:k:(length(m)-rem(length(m),k))
        [txNo, OK] = deal(0);
        frame = m(j:j+k-1);                                 %formação dos blocos
        
        % Rayleigh não-codificado
        bpsk = 2*frame-1;                                   %bpsk
        noise = sqrt(0.5*N0)*(randn(1,k)+1i*randn(1,k));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,k); %rayleigh canal semi-estático
        y = bpsk.*h + noise;                                %transmitido
        err = sum((y./h>0)~= frame) + err;                  %soma de bits errados
        
        while(~OK && txNo < 3)                
            % Codificação Hamming (15,11)
            frame_enc = encode(frame,n,k,'hamming/binary');     %codificação
            bpsk_enc = 2*frame_enc-1;                           %bpsk
            noise = sqrt(0.5*N0)*(randn(1,n)+1i*randn(1,n));	%awgn independe do canal
            h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,n); %rayleigh canal semi-estático
            y = bpsk_enc.*h + noise;                            %sinal transmitido
            
            y_dem = y./h>0;                                     %demodulado
            y_dec = decode(y_dem,n,k,'hamming/binary');         %decodificado
            
            err_enc(txNo+1) = sum(y_dec ~= frame) + err_enc(txNo+1);    %erros acumulados

            if(sum(y_dec ~= frame) > 0)
                txNo = txNo + 1;
            else
                OK = 1;
            end
        end
        
        nBits_enc = nBits_enc + 1;
    end
    
    ber(i) = err/nBits;
    
    for j=1:3
        ber_enc(j,i) = err_enc(j)/(nBits_enc*k);
    end
end

fer = 1-(1-ber).^n;
fer_enc = 1-(1-ber_enc).^n;

plot_slow = figure('NumberTitle', 'off','units','normalized','outerposition',[0 0 1 1],'Name', 'Slow Fading');
semilogy(EbN0_db,fer_enc(1,:),'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,fer_enc(2,:),'-o','Linewidth', 2);
semilogy(EbN0_db,fer_enc(3,:),'-o','Linewidth', 2);
semilogy(EbN0_db,fer,'-d','Linewidth', 2);
legend('Hamming (15,11) - 0 Retransmissões', 'Hamming (15,11) - 1 Retransmissão', ...
    'Hamming (15,11) - 2 Retransmissões', 'Não-codificado','location', 'best');
title('Desempenho BPSK com codificação Hamming (15,11) e Retransmissão HARQ Tipo-I Simples em Canal Lento com Desvanecimento Rayleigh');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0)-1 1e3/nBits 1]);

print(plot_slow,'t6_q4-Thais','-dpng')

close(bar);

%% Q5
%{
Segundo Goldsmith, as técnicas de combinação de pacotes mais comuns são:
Selection Combining, Switch and Stay Combining e Maximal Ratio Combining.

Selection Combining: o combinador retorna o sinal com a maior SNR;

Switch and Stay Combining: o combinador escolhe um canal no qual o sinal
recebido possua uma SNR acima de um certo limite e se mantém nele até que o
sinal possua uma SNR menor, de maneira que o combinador escolherá outro
canal a partir de um critério específico para receber os pacotes.

Maximal Ratio Combining: a combinação é dada pela soma de todos os canais
com um determinado peso. Esta técnica também é sugerida por David Chase
(1985) no artigo entitulado "Code Combining - A Maximum-Likelihood Decoding
Approach for Combining an Arbitrary Number of Noisy Packets", no qual o
autor destaca dois fatores cruciais para a eficiência da técnica de
combinação de pacotes: a presença de um decodificador de
máxima-verossimilhança e a utilização de pesos em cada pacote recebido para
se estimar a sua confiabilidade.
%}

%% Q6
%{
A AMC consiste em técnicas que visam manter uma Eb/N0 constante variando a potência de transmissão, taxa de transmissão, tamanho da constelação, taxa de codificação e suas combinações. Assim, tais esquemas de AMC fornecem uma alta eficiência espectral por transmitir em altas taxas em condições favoráveis do canal, sendo seu desempenho melhorado ao ser combinado com técnicas de diversidade espacial.
%}
