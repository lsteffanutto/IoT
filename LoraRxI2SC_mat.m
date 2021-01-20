function [bitestLoRa, symbolesEstLoRa]=LoraRxI2SC_mat(sigRx,PreambleLength,alpha,twoSF,Ns,chirp_brut,freq_axis)
    t0 = PreambleLength + 1;
    N = twoSF; % Points FFT

    sigRx  = sigRx(t0:end).';
    Mat_sigRx = reshape(sigRx,N,Ns);
    Mat_chirp_brut = repmat(chirp_brut',1,Ns);

    S = Mat_sigRx.*Mat_chirp_brut;
    S_FFT = abs(fft(S,N)).^2;

    [~,symbolesEstLoRa] = max(S_FFT);
    symbolesEstLoRa = (N+1 - symbolesEstLoRa)';
    symbolesEstLoRa(symbolesEstLoRa==N) = 0;

    % BITS
    bitestLoRa = bin2gray(symbolesEstLoRa,'pam',N);
    bitestLoRa = de2bi(bitestLoRa,log2(N),'left-msb');
    bitestLoRa = bitestLoRa>0.5;
end