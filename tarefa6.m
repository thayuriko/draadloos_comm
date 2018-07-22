%% Q2 - BER e FER do BPSK com c�digos com desvanecimento Rayleigh (canais r�pido e lento)

%---------------FAST FADING---------------%
clear all; close all; clc;

rand('state',0); randn('state',0);
bar = waitbar(0,'Setting up...');

EbN0_db = 1:30;
nBits = 10^6;
EbN0 = 10.^(EbN0_db/10);
ber_awgn = zeros(1,length(EbN0));
n = 15; k = 11;

disp('Desempenho C�digo de Hamming (15,11) p/ BPSK');

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait (fast fading)');
    N0 = 1/EbN0(i);
    m = rand(1,nBits)>0.5;  %gera��o de bits aleat�rios
    [nBits_enc, err_enc] = deal(0);
    bpsk = 2*m-1;           %modula��o bpsk sem codifica��o
    
    noise = sqrt(0.5*N0)*(randn(1,nBits)+1i*randn(1,nBits)); %awgn
    h = sqrt(0.5)*(randn(1,nBits)+1i*randn(1,nBits)); %rayleigh
    
    y = bpsk.*h + noise; %sinal transmitido
    
    err = sum((y./h>0)~= m);
    
    for j=1:k:(length(m)-rem(length(m),k))
        frame = m(j:j+k-1); %forma��o dos blocos
        frame_enc = encode(frame,n,k,'hamming/binary'); %codifica��o
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
legend('Codificado - Hamming (15,11)', 'N�o-codificado','location', 'best');
ylabel('BER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

subplot(1,2,2);
semilogy(EbN0_db,fer_enc,'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,fer,'-d','Linewidth', 2);
legend('Codificado - Hamming (15,11)', 'N�o-codificado','location', 'best');
%title('FER x EbN0');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);
close(bar);

a = axes;
t1 = title('Performance BPSK com desvanecimento Rayleigh em Canal R�pido','FontSize',16);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');


%---------------SLOW FADING---------------%
bar = waitbar(0,'Setting up...');
ber_awgn = zeros(1,length(EbN0));

for i=1:length(EbN0)
    waitbar(EbN0_db(i)/length(EbN0),bar,'Please wait (slow fading)');
    N0 =  1/EbN0(i);
    m = rand(1,nBits)>0.5;  %gera��o de bits aleat�rios
    [nBits_enc, err_enc, err] = deal(0);
    
    for j=1:k:(length(m)-rem(length(m),k))
        frame = m(j:j+k-1);                                 %forma��o dos blocos
        
        % Rayleigh n�o-codificado
        bpsk = 2*frame-1;                                   %bpsk
        noise = sqrt(0.5*N0)*(randn(1,k)+1i*randn(1,k));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,k); %rayleigh canal semi-est�tico
        y = bpsk.*h + noise;                                %transmitido
        err = sum((y./h>0)~= frame) + err;                  %soma de bits errados
        
        % Codifica��o Hamming (15,11)
        frame_enc = encode(frame,n,k,'hamming/binary');     %codifica��o
        bpsk_enc = 2*frame_enc-1;                           %bpsk
        noise = sqrt(0.5*N0)*(randn(1,n)+1i*randn(1,n));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,n); %rayleigh canal semi-est�tico
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
legend('Codificado - Hamming (15,11)', 'N�o-codificado','location', 'best');
ylabel('BER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

subplot(1,2,2);
semilogy(EbN0_db,fer_enc,'-o','Linewidth', 2); hold on;
semilogy(EbN0_db,fer,'-d','Linewidth', 2);
legend('Codificado - Hamming (15,11)', 'N�o-codificado','location', 'best');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0) 1e2/nBits 1]);

a = axes;
t1 = title('Performance BPSK com desvanecimento Rayleigh em Canal Lento','FontSize',16);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');

print(plot_slow,'t6_q2_slow-Thais','-dpng')
print(plot_fast,'t6_q2_fast-Thais','-dpng')

close(bar);
%% Q4 - Simular a FER do BPSK com c�digos considerando 0, 1 e 2 retransmiss�es (HARQ Tipo-I Simples)
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
    m = rand(1,nBits)>0.5;  %gera��o de bits aleat�rios
    [err, nBits_enc] = deal(0);
    err_enc = zeros(1,3);
    
    for j=1:k:(length(m)-rem(length(m),k))
        [txNo, OK] = deal(0);
        frame = m(j:j+k-1);                                 %forma��o dos blocos
        
        % Rayleigh n�o-codificado
        bpsk = 2*frame-1;                                   %bpsk
        noise = sqrt(0.5*N0)*(randn(1,k)+1i*randn(1,k));	%awgn independe do canal
        h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,k); %rayleigh canal semi-est�tico
        y = bpsk.*h + noise;                                %transmitido
        err = sum((y./h>0)~= frame) + err;                  %soma de bits errados
        
        while(~OK && txNo < 3)                
            % Codifica��o Hamming (15,11)
            frame_enc = encode(frame,n,k,'hamming/binary');     %codifica��o
            bpsk_enc = 2*frame_enc-1;                           %bpsk
            noise = sqrt(0.5*N0)*(randn(1,n)+1i*randn(1,n));	%awgn independe do canal
            h = sqrt(0.5)*repmat(randn(1,1)+1i*randn(1,1),1,n); %rayleigh canal semi-est�tico
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
legend('Hamming (15,11) - 0 Retransmiss�es', 'Hamming (15,11) - 1 Retransmiss�o', ...
    'Hamming (15,11) - 2 Retransmiss�es', 'N�o-codificado','location', 'best');
title('Desempenho BPSK com codifica��o Hamming (15,11) e Retransmiss�o HARQ Tipo-I Simples em Canal Lento com Desvanecimento Rayleigh');
ylabel('FER'); xlabel('EbN0 (dB)');
grid on; axis([1 length(EbN0)-1 1e3/nBits 1]);

print(plot_slow,'t6_q4-Thais','-dpng')

close(bar);

%% Q5
%{
Segundo Goldsmith, as t�cnicas de combina��o de pacotes mais comuns s�o:
Selection Combining, Switch and Stay Combining e Maximal Ratio Combining.

Selection Combining: o combinador retorna o sinal com a maior SNR;

Switch and Stay Combining: o combinador escolhe um canal no qual o sinal
recebido possua uma SNR acima de um certo limite e se mant�m nele at� que o
sinal possua uma SNR menor, de maneira que o combinador escolher� outro
canal a partir de um crit�rio espec�fico para receber os pacotes.

Maximal Ratio Combining: a combina��o � dada pela soma de todos os canais
com um determinado peso. Esta t�cnica tamb�m � sugerida por David Chase
(1985) no artigo entitulado "Code Combining - A Maximum-Likelihood Decoding
Approach for Combining an Arbitrary Number of Noisy Packets", no qual o
autor destaca dois fatores cruciais para a efici�ncia da t�cnica de
combina��o de pacotes: a presen�a de um decodificador de
m�xima-verossimilhan�a e a utiliza��o de pesos em cada pacote recebido para
se estimar a sua confiabilidade.
%}

%% Q6
%{
A AMC consiste em t�cnicas que visam manter uma Eb/N0 constante variando a pot�ncia de transmiss�o, taxa de transmiss�o, tamanho da constela��o, taxa de codifica��o e suas combina��es. Assim, tais esquemas de AMC fornecem uma alta efici�ncia espectral por transmitir em altas taxas em condi��es favor�veis do canal, sendo seu desempenho melhorado ao ser combinado com t�cnicas de diversidade espacial.
%}
