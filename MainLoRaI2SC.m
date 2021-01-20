clear ;
%close all;

%% Programme principal Tx/Rx LoRa

%% Param?tres d'initialisations

SF = 7;            % Spreading Factor = nb de bits/symboles
M = 2^(SF);
BwL = 125e3;        % Largeur de bande du signal (bande balayee par le chirp)
Ts = M/BwL;         % definition LoRa
Db = SF/Ts;         % debit binaire
P = 1;              % Puissance en watt du signal transmis

alpha = 1;          % facteur de surechantillonnage des chirps (pour pouvoir simuler des desynchro temps)

Te = Ts / (alpha*M); % Periode d'echantillonnage
Fe = 1/Te;

NbOctet=255;            % Nombre d'octect dans la payload
Ns=floor(NbOctet*8/SF); % Nombre de symbole dans la payload

NbEch = Ts*Fe;          % Nombre d'echantillons dans la payload
NbPaquet=100;            % Nombre de paquet transmis (utiliser pour tracer des courbes de performances BER par exemple)
SNRdB = 30;%-30:1:-10;             % Rapport signal sur bruit au recepteur defini comme le rapport entre la puissance du signal recu divis? par la puissance du bruit
NF = -8;                % Facteur de bruit typique d'un recepteur 
Sensitivity = -174 + 10*log10(BwL) + SNRdB + NF;

BER = zeros(1,length(SNRdB));
SER = zeros(1,length(SNRdB));
PER = zeros(1,length(SNRdB));

%% Iteration Paquets
for numSNR = 1 : length(SNRdB)
    fprintf('Iteration %d/%d, SF = %d, SNR = %ddB \n',numSNR,length(SNRdB),SF, SNRdB(numSNR))
    for numPaquet = 1:NbPaquet
        
        %% Emetteur LoRa
        [sig, bits, PreambleLength, chirp_brut, freq_axis, symboles]=LoraTxI2SC(SF,BwL,alpha,Ns,M,Ts,Te);
        
        %% Canal
        h = 1; % suppose non selectif en frequence
        sigRx = filter(h,1,sig);
        
        %% Recepteur LoRa
        % On considere que les signaux recus sont synchronises temporellement
        % et frequentiellement
        
        Psig = mean(abs(sigRx).^2); % puissance du signal recu
        sigRx = sqrt(P/Psig)*sigRx; % On ajuste la puissance du signal recu pour quelle soit egale a P
        
        Pb = P/10.^(SNRdB(numSNR)/10);   % puissance du bruit
        
        bruit = sqrt(Pb/2)*(randn(size(sigRx))+1i*randn(size(sigRx))); % Generation du bruit
        
        % Spectrogramme
%         figure
%         spectrogram(sigRx,32,16,16,Fe,'yaxis')
%         colormap('jet');
%         title('Spectrogramme d un signal non bruite avec SF = 7')
        
        sigRx = sigRx + bruit;
        
        [bitestLoRa, symbolesEstLoRa]=LoraRxI2SC(sigRx,PreambleLength,alpha,2^SF,Ns,chirp_brut,freq_axis);
        
        
        
        %% Calcul BER, SER et PER LoRa
        NbBitsFaux = sum(abs(reshape(bits.',1,[])-bitestLoRa));
        NbBitsTotal = length(bitestLoRa);
        BER(numSNR) = BER(numSNR) + NbBitsFaux/NbBitsTotal;
        SER(numSNR) = SER(numSNR) + mean(symbolesEstLoRa~=symboles.');        
        if NbBitsFaux>0
            PER(numSNR)=PER(numSNR)+1;
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

save('SF_12_NewLoRa','BER','SER','PER','SNRdB');



