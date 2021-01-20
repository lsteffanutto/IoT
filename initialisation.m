function [M,T,Db,Ts,Fe,Ns,NbEch] = initialisation(SF,alpha,BwL)

    M = 2^(SF);
    T = M/BwL;         % definition LoRa
    Db = SF/T;         % debit binaire

    Ts = T / (alpha*M); % Periode d'echantillonnage
    Fe = 1/Ts;

    NbOctet=255;            % Nombre d'octect dans la payload
    Ns=floor(NbOctet*8/SF); % Nombre de symbole dans la payload, nombre de bit par symboles SF

    NbEch = T*Fe;          % Nombre d'echantillons dans la payload
end

