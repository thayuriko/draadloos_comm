%% Q.1
%{
Sabendo-se que a probabilidade de erro de bit da modula��o BPSK � dada por:

$P_{b}=Q\begin{pmatrix}
\frac{^{d_{min}}}{\sqrt{2N_{0}}}
\end{pmatrix}$

Sendo que $d_{min}$ corresponde � dist�ncia entre bits e � igual a 2A, pode-se relacionar A � energia de bit por:

$E_{b}=\int_{0}^{T_{b}}s_{1}^{2}(t)dt=\int_{0}^{T_{b}}s_{2}^{2}(t)dt=\int_{0}^{T_{b}}A^{2}g^{2}(t)cos^{2}(2\pi f_{c}t)dt=A^{2}$

Ent�o a m�nima dist�ncia entre bits corresponde � $2A = 2\sqrt{E_{b}}$

Portanto:

$P_{b}=Q\begin{pmatrix}
\frac{^{2\sqrt{E_{b}}}}{\sqrt{2N_{0}}}
\end{pmatrix}=Q\begin{pmatrix}
\sqrt{\frac{^{2E_{b}}}{N_{0}}}
\end{pmatrix}=Q(\sqrt{2\gamma_{b}})$

Em um canal com desvanecimento do tipo Rayleigh calcula-se a probabilidade m�dia de erro do sinal, utilizada quando $T_{c}\approx T_{s}$. Assim, como a SNR � praticamente constante no tempo, a probabilidade m�dia de erro � calculada pela integral da probabilidade de erro do sinal em um canal AWGN com uma distribui��o Rayleigh, a qual � dada por:

$p_{\gamma_{b}}(\gamma)=\frac{1}{\overline{\gamma_{b}}}e^{-\frac{\gamma}{\overline{\gamma_{b}}}}$

A qual, integrada com a probabilidade de erro de bit para a modula��o BPSK, resulta na probabilidade m�dia de erro de bit para o desvanecimento Rayleigh:

$\overline{P_{\gamma}}=\frac{1}{2}\begin{bmatrix}
1-\sqrt\frac{\overline{\gamma_{b}}}{1+\overline{\gamma_{b}}}
\end{bmatrix}\approx \frac{1}{4\overline{\gamma_{b}}}$
%}

%% Q.2
%{
A partir a equa��o (6.61), a express�o geral da probabilidade m�dia de erro de bit em um canal com desvanecimento Rayleigh pode ser aproximada como:

$\overline{P_{b}}\approx \frac{\alpha_{M}}{2} \left (
1-\sqrt{\frac{0.5\beta_{M} \overline{\gamma_{s}}}{1+0.5\beta_{M}\overline{\gamma_{b}}}}
\right )$

sendo $\alpha_{M}$ e $\beta_{M}$ tais quais que $P_{b}\approx \alpha_{M} Q(\sqrt{\beta_{M} \gamma_{b}})$, variando para cada tipo de modula��o.

Portanto, como a express�o geral da probabilidade de erro de bit para a M-PSK � dada por:

$P_{b}\approx \frac{2}{\log_{2}M}Q(\sqrt{2\gamma_{b}\log_{2}M}sin(\pi /M))$ (6.18)

Temos que $\alpha_{M} = \frac{2}{\log_{2}M}$ e $\beta_{M} = 2\log_{2}M\sin^2(\pi /M))$. Logo, a probabilidade de erro de bit para M-PSK com desvanecimento Rayleigh se torna:

$\overline{P_{b}}\approx \frac{1}{\log_{2}M} \left (
1-\sqrt{\left (\frac{\log_{2}M.\overline{\gamma_{b}}}{1+\log_{2}M.\sin^2(\pi /M)\overline{\gamma_{b}}} \right )}sin(\pi /M)
\right )$


Para o caso da modula��o M-QAM retangular, sua probabilidade de erro de bit � dada por:

$P_{b}\approx \frac{2(\sqrt{M}-1)}{\sqrt{M}\log_{2}M}Q
\left(
\sqrt{\frac{3\overline{\gamma_{b}}\log_{2}M}{M-1}}
\right)$, (6.24)

sendo os par�metros $\alpha_{M} = \frac{2(\sqrt{M}-1)}{\sqrt{M}\log_{2}M}$ e $\beta_{M} = \frac{3\log_{2}M}{M-1}$.

Assim, obt�m-se a $P_{s}$ m�dia para o desvanecimento Rayleigh e modula��o M-QAM:

$\overline{P_{b}}\approx \frac{(\sqrt{M}-1)}{\sqrt{M}\log_{2}M}
\left(
1-\sqrt{\frac{(M-1)1.5\overline{\gamma_{b}}\log_{2}M}{(M-1)^2+1.5\overline{\gamma_{b}}\log_{2}M}}
\right)$
%}

%% Q.3
clear all; close all; clc;

Pt = 4e-3; Gt_db = 5; Gr_db = 8; Pl_db = 100; margem_db = 5; 
N0_db = -204; BER = 1e-3; B = 5e6; alpha = 0.25; Bmax = 0;

Pt_db = 10*log10(Pt);
Pl = 10^(Pl_db/10);
N0 = 10^(N0_db/10);
Pr_db = Pt_db + (Gt_db + Gr_db) - Pl_db - margem_db;
Pr = 10^(Pr_db/10);
M = 2;

EbN0 = ((qfuncinv(BER*log2(M)/2)/sin(pi/M))^2)/(2*log2(M));
EbN0_db = 10*log10(EbN0);

EbN0_h = ((1-BER*2)^2)/(1 - (1-BER*2)^2);
EbN0_db_h = 10*log10(EbN0_h);

%{
Para os par�metros estipulados pelo exerc�cio, a nova SNR calculada para um
canal com desvanecimento Rayleigh seria de 24dB, sendo quase 17dB maior do
que a SNR necess�ria para obter a BER exigida no problema em quest�o. Para
alterar a SNR, seria preciso alterar a ordem da modula��o e a BER
requerida, portanto, para este cen�rio especificado, n�o seria poss�vel
transmitir pacotes em um canal com desvanecimento Rayleigh.
%}

%% Q.4
clear all; close all; clc;

rand('state',0);
randn('state',0);

EbN0_db = 0:40;  % vetor de rela��o sinal-ru�do (SNR) em dB
nBits = 10^7;    % quantidade de bits transmitidos
Eb = 1;          % energia de bit   
EbN0 = 10.^(EbN0_db/10);   % convers�o da SNR de dB para linear;
ber_awgn = zeros(1,length(EbN0)); % inicializa��o do vetor de BER
m_naka = 3;

for i=1:length(EbN0)
    N0 =  Eb/EbN0(i);
    m = rand(1,nBits)>0.5;
    x = 2*m-1;

    n = sqrt(0.5*N0)*(randn(1,nBits)+1i*randn(1,nBits));
    h = sqrt(0.5)*(randn(1,nBits)+1i*randn(1,nBits)); %gaussiana na amplitude e fase (divide por sqrt(1/2) pra ter amplitude unit�ria)
    r = nak_m(m_naka,1,nBits); %randraw
    
    y0 = x + n;
    y1 = x.*r + n;
    y2 = x.*h + n;
    
    err0 = sum((y0>0)~= m);
    err1 = sum((y1./r>0)~= m);
    err2 = sum((y2./h>0)~= m);
    
    ber_awgn(i) = err0/nBits;
    ber_naka(i) = err1/nBits;
    ber_ray(i) = err2/nBits;
    
    Pb_ray(i) = 1/pi*integral(@(phi)mgf('rayl',phi,m_naka,EbN0(i)),0,pi/2);
    Pb_naka(i) = 1/pi*integral(@(phi)mgf('naka',phi,m_naka,EbN0(i)),0,pi/2);
    Pb_awgn(i) = qfunc(sqrt(2*EbN0(i)));
end

semilogy(EbN0_db,Pb_ray,'Linewidth', 2); hold on;
semilogy(EbN0_db,ber_ray,'o','Linewidth', 2);
semilogy(EbN0_db,Pb_naka,'Linewidth', 2);
semilogy(EbN0_db,ber_naka,'o','Linewidth', 2);
semilogy(EbN0_db, Pb_awgn,'Linewidth', 2);
semilogy(EbN0_db, ber_awgn, 'o','Linewidth', 2);

legend('Pb Rayleigh', 'BER Rayleigh', 'Pb Nakagami m = 3', 'BER Nakagami m = 3', 'Pb AWGN', 'BER AWGN');
title('Desempenho do canal para as distribui��es AWGN, Rayleigh e Nakagami-m');
ylabel('Pb ou BER'); xlabel('EbN0 (dB)'); axis([0 40 100/nBits 1]);
grid on;