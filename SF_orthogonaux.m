clear ; clc
close all;
clc

hold all
%% Programme principal Tx/Rx LoRa

%% Parametres
SF1 = 7;
SF2 = 8;
alpha = 1;          % facteur de surechantillonnage des chirps (pour pouvoir simuler des desynchro temps)
BwL = 125e3;           % Largeur de bande du signal (bande balayee par le chirp)

[M1,T1,Db1,Ts1,Fe1,Ns1,NbEch1] = initialisation(SF1,alpha,BwL);
[M2,T2,Db2,Ts2,Fe2,Ns2,NbEch2] = initialisation(SF2,alpha,BwL);

P = 14;              % Puissance en watt du signal transmis
NbOctet=255;            % Nombre d'octect dans la payload
NbPaquet=1e1;            % Nombre de paquet transmis (utiliser pour tracer des courbes de performances BER par exemple)
SNRdB = -20:1:5;          % Rapport signal sur bruit au recepteur defini comme le rapport entre la puissance du signal recu divis??? par la puissance du bruit

NF = -8;                  % Facteur de bruit typique d'un recepteur
Sensitivity = -174 + 10*log10(BwL) + SNRdB + NF;


%% Orthogonalite
PR1 = -10:2:10;
lPR1 = length(PR1);


%% Initialisation des tableaux
BER1 = zeros(lPR1 ,length(SNRdB));   % bit
BER2 = zeros(lPR1 ,length(SNRdB));   % bit

figure(1)
%figure(2)

for PR=1:lPR1
    for numSNR = 1 : length(SNRdB)
        fprintf('Iteration %d/%d, SF = %d, SNR = %ddB \n',numSNR,length(SNRdB),SF1, SNRdB(numSNR))
        for numPaquet = 1:NbPaquet
            
            %% Emetteur LoRa
            [Ssf1, bits1, PreambleLength1, chirp_brut1, ~, symboles1]         = LoraTxI2SC(SF1,BwL,alpha,Ns1,M1,T1,Ts1);
            [Ssf2, bits2, PreambleLength2, chirp_brut2, freq_axis, symboles2] = LoraTxI2SC(SF2,BwL,alpha,Ns2,M2,T2,Ts2);
            
            Psf1 = 1;
            Psf2 = Psf1.*(10^-(PR1(PR)/10));

            Psf2_real = mean(abs(Ssf2).^2);
            Psf1_real = mean(abs(Ssf1).^2);
            sigRx = [sqrt(Psf1/Psf1_real).*Ssf1 + sqrt(Psf2/Psf2_real).*Ssf2(1:length(Ssf1)) sqrt(Psf2/Psf2_real).*Ssf2(length(Ssf1)+1:end)];
            
            %% Recepteur LoRa
            % On considere que les signaux recus sont synchronises temporellement et frequentiellement
            
            Psig = mean(abs(sigRx).^2); % puissance du signal 1 recu
            sigRx = sqrt(P/Psig)*sigRx; % On ajuste la puissance du signal recu pour quelle soit egale a P
            
            Pb = P/10.^(SNRdB(numSNR)/10);   % puissance du bruit
            
            bruit = sqrt(Pb/2)*(randn(size(sigRx))+1i*randn(size(sigRx))); % Generation du bruit
            
            sigRx = sigRx + bruit;
            
            %% Notre decodeur LoRa
            
          [bitestLoRa1, symbolesEstLoRa1]=LoraRxI2SC_mat(sigRx(1:length(Ssf1)),PreambleLength1,alpha,2^SF1,Ns1,chirp_brut1,freq_axis);
            
          [bitestLoRa2, symbolesEstLoRa2]=LoraRxI2SC_mat(sigRx,PreambleLength2,alpha,2^SF2,Ns2,chirp_brut2,freq_axis);
            
           %% Calcul BER, SER et PER LoRa
            
            NbBitsFaux = sum(abs(reshape(bits1.',1,[])-reshape(bitestLoRa1.',1,[]) ));
            NbBitsTotal = length(bitestLoRa1);
            BER1(PR,numSNR) = BER1(PR,numSNR) + NbBitsFaux/NbBitsTotal;
            
            
            NbBitsFaux = sum(abs(reshape(bits2.',1,[])-reshape(bitestLoRa2.',1,[])));
            NbBitsTotal = length(bitestLoRa2);
            BER2(PR,numSNR) = BER2(PR,numSNR) + NbBitsFaux/NbBitsTotal;
            
                       
        end
    end
    
    BER1 = BER1 / NbPaquet;
    BER2 = BER2 / NbPaquet;

    
    color = {''};
    
    figure(1)
    semilogy(Sensitivity,[BER1(PR,:)]);
    xlabel('Sensitivity (dBm)');
    ylabel('BER1');
    legend('BER');
    title(sprintf('Performane LoRa communication SF1 = %d',SF1))
    grid on;
    hold on;
    
%     figure(2)
%     semilogy(Sensitivity,[BER(PR,:)]);
%     hold on;
%     semilogy(Sensitivity,[BER2(PR,:)]);
%     xlabel('Sensitivity (dBm)');
%     ylabel('BER2');
%     legend('BER');
%     title(sprintf('Performane LoRa communication SF2 = %d',SF2))
%     grid on;
%     hold on;
    
    
end
legend('BER ref','PR1 = -10','PR1 = -8','PR1 = -6','PR1 = -4','PR1 = -2','PR1 = 0','PR1= 2','PR1 = 4','PR1 = 6','PR1 = 8','PR1 = 10');
%% Sauvegarde
save(['BER SF1=',num2str(SF1),'SF2=',num2str(SF2),'orthogonalite_SF'],'Sensitivity','PR1','BER1','BER2','SNRdB');





