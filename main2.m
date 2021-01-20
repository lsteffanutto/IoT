clear;
close all;
clc
%% Programme principal Tx/Rx LoRa

%% Parametres d'initialisations

SF = 7;
SNRdB = -30:1:-10;             
BER = zeros(6,length(SNRdB));
SER = zeros(6,length(SNRdB));
PER = zeros(6,length(SNRdB));
alpha = 1;          % facteur de surechantillonnage des chirps (pour pouvoir simuler des desynchro temps)
P = 1;              % Puissance en watt du signal transmis
NF = -8;                % Facteur de bruit typique d'un recepteur
for SF = 7%:12 % Spreading Factor = nb de bits/symboles
    
    M = 2^(SF);
    BwL = 125e3;        % Largeur de bande du signal (bande balayee par le chirp)
    T = M/BwL;          % definition LoRa
    Db = SF/T;          % debit binaire    
    
    Ts = T / (alpha*M); % Periode d'echantillonnage
    Fe = 1/Ts;
    
    NbOctet=255;            % Nombre d'octect dans la payload
    Ns=floor(NbOctet*8/SF); % Nombre de symbole dans la payload
    
    NbEch = T*Fe;          % Nombre d'echantillons dans la payload
    NbPaquet=1e2;            % Nombre de paquet transmis (utiliser pour tracer des courbes de performances BER par exemple)
    
    Sensitivity = -174 + 10*log10(BwL) + SNRdB + NF;
    
    %% Iteration Paquets
    
    for numSNR = 1 : length(SNRdB)
        fprintf('Iteration %d/%d, SF = %d, SNR = %ddB \n',numSNR,length(SNRdB),SF, SNRdB(numSNR))
        for numPaquet = 1:NbPaquet
            
            %% Emetteur LoRa
            [sig, bits, PreambleLength, chirp_brut, freq_axis, symboles]=LoraTxI2SC(SF,BwL,alpha,Ns,M,T,Ts);
            
            %% Canal
            h = 1; % suppose non selectif en frequence
            sigRx = filter(h,1,sig);
            
            %% Recepteur LoRa
            % On considere que les signaux recus sont synchronises temporellement
            % et frequentiellement
            
            Psig = mean(abs(sigRx).^2); % puissance du signal reeu
            sigRx = sqrt(P/Psig)*sigRx; % On ajuste la puissance du signal recu pour quelle soit egale a P
            
            Pb = P/10.^(SNRdB(numSNR)/10);   % puissance du bruit
            
            bruit = sqrt(Pb/2)*(randn(size(sigRx))+1i*randn(size(sigRx))); % Generation du bruit
            
            sigRx = sigRx + bruit;
            
            [bitestLoRa, symbolesEstLoRa]=LoraRxI2SC(sigRx,PreambleLength,alpha,2^SF,Ns,chirp_brut,freq_axis);
            
                      
            %% Calcul BER, SER et PER LoRa
            NbBitsFaux = sum(abs(reshape(bits.',1,[])-bitestLoRa));
            NbBitsTotal = length(bitestLoRa);
            BER(SF-6,numSNR) = BER(SF-6,numSNR) + NbBitsFaux/NbBitsTotal;
            SER(SF-6,numSNR) = SER(SF-6,numSNR) + mean(symbolesEstLoRa~=symboles.');
            if NbBitsFaux>0
                PER(SF-6,numSNR)=PER(SF-6,numSNR)+1;
            end
            
        end
    end
    
    
    BER = BER / NbPaquet;
    SER = SER / NbPaquet;
    PER = PER / NbPaquet;

    figure(1)
    semilogy(Sensitivity,[BER;SER;PER]);
    xlabel('Sensitivity (dBm)');
    ylabel('BER / SER / PER');
    legend('BER','SER','PER');
    title(sprintf('Performane LoRa communication SF = %d',SF))
    grid on;
    hold on;


    
end

save('Parti_main2.mat','BER','SER','PER')
 

%% spectrogrammes

[sig, bits, PreambleLength, chirp_brut, freq_axis, symboles]=LoraTxI2SC(SF,BwL,alpha,Ns,M,T,Ts);

% % Spectrogramme
% spectrogram(sig,32,16,16,Fe,'yaxis');
% title('Spectrogramme du signal emis bruit avec ZF=7')
% 
% % Spectrogramme
% spectrogram(sigRx,32,16,16,Fe,'yaxis');
% title('Spectrogramme du signal recu sans avec bruit ZF=7')



