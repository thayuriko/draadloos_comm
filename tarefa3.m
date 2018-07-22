%% Q1
clear all; close all; clc;

pr = 1; j = 1; h = 0:0.1:2.5; K = [0 2 5 10]; m = [0.5 1 2 5];

phRice = zeros(1,length(K));
Markers = {'d','o','*','^'};
Marker_Counter = 1;

legendCell = cell(length(K), 2);

for i=1:length(K)
    z = 2.*h.*sqrt(K(i)*(K(i)+1)/pr);
    phRice =  2.*h*(K(i)+1).*exp(-K(i)-(K(i)+1).*h.^2/pr).*besseli(0,z)/pr;
    
    phNaka = (2*m(i)^m(i).*h.^(2*m(i)-1)).*exp(-m(i).*h.^2/pr)/gamma(m(i))*pr.^m(i);
    
    figure(1)
    plot(phRice,strcat('-',Markers{Marker_Counter}));
    hold on;
    
    figure(2)
    plot(phNaka,strcat('-',Markers{Marker_Counter}));
    hold on;    
    
    Marker_Counter = Marker_Counter+1;
    
    if m(i) == 1
        legendCell{i,2} = ['\mu = 1 (Rayleigh)'];
    else
        legendCell{i,2} = ['\mu = ' num2str(m(i))];
    end
    
    if K(i) == 0
        legendCell{i,1} = ['K = 0 (Rayleigh)'];
    else
        legendCell{i,1} = ['K = ' num2str(K(i))];
    end
end

figure(1)
legend(legendCell{:,1});
title('PDF Rice p/ diferentes valores de K');
ylabel('p_h(h)'); xlabel('h');
grid on;

figure(2)
legend(legendCell{:,2});
title('PDF Nakagami p/ diferentes valores de \mu');
ylabel('p_h(h)'); xlabel('h');
grid on;

figure(3)
z = 0;
phRay =  2.*h.*exp(-h.^2/pr).*besseli(0,z)/pr;
plot(phRay);
legend('Rayleigh');
title('PDF Rayleigh');
ylabel('p_h(h)'); xlabel('h');
grid on;
%rayleigh é uma gaussiana real + uma gaussiana imaginaria
%abs(randn(1,1e6) + 1i*randn(1,1e6))
%

%% Q3
% existem funções prontas no matlab para isso:
% randn: gaussiana
% rand: uniforme
% é necessário ter uma amostra da distribuição uniforme, sua CDF: inverse
% transform sampling



%% Q4
clear all; close all; clc;

v = 120/3.6; fc = 900e6; B = 200e3; sigma = 25e-6;
c = 3e8;

lambda = c/fc;
% parametros de dispersão de frequencia: espalhamento doppler

Bd = v/lambda; %espalhamento doppler
% compara-se com o tempo de coerência
Tc = 0.423/Bd; %em segundos
%como o tempo é de 4ms, a duração de um time slot é 0.5ms, logo como o Tc >
%Ts, a seletividade é lenta

Bc = 1/(5*sigma);
%como a banda é de 200k e Bc =  8k, o canal muda mto rápido na frequência,
%logo é seletivo na frequência
fprintf('Q4. ');

%% Q5ans
clear all; close all; clc;

atraso = [0 .3 3.5 4.4 9.5 12.7]; %atraso/potencia
potencia_db = [0 -12 -4 -7 -15 -22];
potencia = 10.^(potencia_db/10);

tau1 = sum(potencia.*atraso)/sum(potencia);
tau2 = sum(potencia.*atraso.^2)/sum(potencia);

sigma = sqrt(tau2 - tau1^2);

Bc = 1/(5*sigma);

fprintf('Q5. A banda de coerência é de %s kHz\n', num2str(Bc*1e3));