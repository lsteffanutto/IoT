clear;
close all;
clc
%% Programme principal Tx/Rx LoRa

%% Parametres d'initialisations
SF1 = 7;
SF2 = 7;
alpha = 1;          % facteur de surechantillonnage des chirps (pour pouvoir simuler des desynchro temps)
P = 1;              % Puissance en watt du signal transmis on est plutot de l'ordre du miliwatt en vrai
BwL = 125e3;           % Largeur de bande du signal (bande balayee par le chirp)

[M1,T1,Db1,Ts1,Fe1,Ns1,NbEch1] = initialisation(SF1,alpha,BwL);
[M2,T2,Db2,Ts2,Fe2,Ns2,NbEch2] = initialisation(SF2,alpha,BwL);

NbPaquet=1e1;            % Nombre de paquet transmis (utiliser pour tracer des courbes de performances BER par exemple)
SNRdB = -20:1:0;             % Rapport signal sur bruit au recepteur defini comme le rapport entre la puissance du signal recu divis??? par la puissance du bruit
SNR = -10;
NF = -8;                % Facteur de bruit typique d'un recepteur 
Sensitivity = -174 + 10*log10(BwL) + SNRdB + NF;


%% Orthogonalite
PR1 = 8;
lPR1 = length(PR1);


%% Initialisation des tableaux
BER1 = zeros(lPR1 ,length(SNRdB));   % bit
SER1 = zeros(lPR1,length(SNRdB));   % symbole
PER1 = zeros(lPR1,length(SNRdB));   % packet

BER2 = zeros(lPR1 ,length(SNRdB));   % bit
SER2 = zeros(lPR1,length(SNRdB));   % symbole
PER2 = zeros(lPR1,length(SNRdB));   % packet

%% Iteration Paquets
figure(1)
figure(2)
for PR=1:lPR1
    for numSNR = 1 : length(SNRdB)
        fprintf('Iteration %d/%d, SF = %d, SNR = %ddB \n',numSNR,length(SNRdB),SF1, SNRdB(numSNR))
        for numPaquet = 1:NbPaquet

            %% Emetteur LoRa
            [Ssf1, ~, PreambleLength1, chirp_brut1, ~, ~]=LoraTxI2SC(SF1,BwL,alpha,Ns1,M1,T1,Ts1);
            [Ssf2, bits2, PreambleLength2, chirp_brut2, freq_axis, symboles2]=LoraTxI2SC(SF2,BwL,alpha,Ns2,M2,T2,Ts2);

            
            Psf1 = 1;
            Psf2 = Psf1.*(10^-(PR1(PR)/10));

            Psf2_real = mean(abs(Ssf2).^2);
            Psf1_real = mean(abs(Ssf1).^2);
            sigRx = sqrt(Psf1/Psf1_real).*Ssf1 + sqrt(Psf2/Psf2_real).*Ssf2;
            
            %% Recepteur LoRa
            % On considere que les signaux recus sont synchronises temporellement
            % et frequentiellement


            Pb = Psf1/10.^(SNRdB(numSNR)/10);   
            bruit = 0;%sqrt(Pb/2)*(randn(size(sigRx))+1i*randn(size(sigRx))); % Generation du bruit

            sigRx = sigRx + bruit ; 

            
            %% Decodage SIC
            % Decodage s1
            [bitestLoRa1, symbolesEstLoRa1]=LoraRxI2SC(sigRx(1:length(Ssf1)),PreambleLength1,alpha,2^SF1,Ns1,chirp_brut1,freq_axis);

            % Reconstruction s1

            [Ssf1_hat, bits1, PreambleLength1, chirp_brut1, freq_axis, symboles1] = LoraTxI2SC_SIC(bitestLoRa1,SF1,BwL,alpha,Ns1,M1,T1,Ts1);
            
            % Mise a la bonne puissance
            Psf1_hat = mean(abs(Ssf1_hat).^2);
            Ssf1_hat = sqrt(Psf1/Psf1_hat)*Ssf1_hat;
            
            % SIC : premiere iteration
            z = sigRx-Ssf1_hat;
            
            test = Ssf1_hat - Ssf1;
            % Decodage s2
            [bitestLoRa2, symbolesEstLoRa2]=LoraRxI2SC(z,PreambleLength2,alpha,2^SF2,Ns2,chirp_brut2,freq_axis);  

                
            %% Calcul BER, SER et PER LoRa
            
            NbBitsFaux = sum(abs(reshape(bits1.',1,[])-bitestLoRa1));
            NbBitsTotal = length(bitestLoRa1);
            BER1(PR,numSNR) = BER1(PR,numSNR) + NbBitsFaux/NbBitsTotal;
            SER1(PR,numSNR) = SER1(PR,numSNR) + mean(symbolesEstLoRa1~=symboles1.');        
            if NbBitsFaux>0
                PER1(PR,numSNR)=PER1(PR,numSNR)+1;
            end
                
          
            NbBitsFaux = sum(abs(reshape(bits2.',1,[])-bitestLoRa2));
            NbBitsTotal = length(bitestLoRa2);
            BER2(PR,numSNR) = BER2(PR,numSNR) + NbBitsFaux/NbBitsTotal;
            SER2(PR,numSNR) = SER2(PR,numSNR) + mean(symbolesEstLoRa2~=symboles2.');        
            if NbBitsFaux>0
                PER2(PR,numSNR)=PER2(PR,numSNR)+1;
            end

        end
    end
    
    BER1 = BER1 / NbPaquet;
    BER2 = BER2 / NbPaquet;
    SER1 = SER1 / NbPaquet;
    SER2 = SER2 / NbPaquet;
    PER1 = PER1 / NbPaquet;
    PER2 = PER2 / NbPaquet;

    color = {''};

    figure(1)
    semilogy(Sensitivity,[BER1(PR,:)]);
    grid on;
    hold on;
    figure(2)
    semilogy(Sensitivity,[BER2(PR,:)]);
    grid on;
    hold on;
    

end
%% End plot

figure(1)
xlabel('Sensitivity (dBm)');
ylabel('BER');
title(sprintf('Performane LoRa communication interferences s1 SF = %d',SF1))

figure(2)
xlabel('Sensitivity (dBm)');
ylabel('BER');
title(sprintf('Performane LoRa communication interferences s2 SF = %d',SF2))


%% Plot
% figure(3)
% hold on
% semilogy(Sensitivity,BER);
% legend([num2str(PR1.');'s1 seul']);
% 
% figure(4)
% hold on
% semilogy(Sensitivity,BER);
% legend([num2str(PR1.');'s2 seul']);



