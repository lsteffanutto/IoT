function [bitestLoRa, symbolesEstLoRa] = LoraRxI2SC(sigRx,PreambleLength,alpha,M,Ns,chirp_brut,freq_axis)
        
    % On enleve le preambule 
    SF = log2(M);
    sig = sigRx(PreambleLength+1:end);
    sig_dechirpe = zeros(1,length(sig));
    N = length(freq_axis);
    BwL = 125e3;        % Largeur de bande du signal (bande balayee par le chirp)
    Ts = M/BwL; 
    Te = Ts / (alpha*M); % Periode d'echantillonnage
    Fe = 1/Te;
    
    % bit symbole
    nb = Ns*SF;
    bitestLoRa = zeros(1,nb);
    symbolesEstLoRa = zeros(1,Ns);

    % Bouclage sur lensemble des symboles
    for i = 1:Ns
        % On prend les N premiers a chaque tour
        sig_N = sig((i-1)*N+1:i*N);
        rm = sig_N.*conj(chirp_brut);
        sig_dechirpe((i-1)*N+1:i*N) = rm;

        % On prend le maximum calcule
        [~, indice] = max(abs(1/sqrt(N)*fft(rm,N)).^2);
        % Symbole estime
        symbolesEstLoRa(i) = mod(M - indice+1, M);

        % passage en bits
        b2g = bin2gray(symbolesEstLoRa(i),'pam',M);  % passage en gray fonction LoraTxI2SC     
        b2g_b = de2bi(b2g,'left-msb');
        len = length(b2g_b);
        if (len ~= SF)
            b2g_b = [zeros(1,SF-len), b2g_b];
        end
        
        bitestLoRa((i-1)*SF+1:i*SF) = b2g_b;            
    end 
    
    %% Affichage du spectrogramme non bruite
% figure
% spectrogram(sig_dechirpe,32,16,16,Fe,'yaxis')
% colormap('jet');
% title('Spectrogramme d un signal dechirpe non bruite avec SF = 7')

end