function [sig, bits, PreambLegnth, chirp_up, freq_axis, symb]=LoraTxI2SC(SF,B,alpha,Ns,M,Ts,Te)

freq_axis=0:M-1;
%% Initialisation
time=[-Ts/2:Te:Ts/2-Te];                % base de temps sur laquelle les chirps sont generes
fc_t=B/(2*Ts)*time;                     % frequence du chirp brut
chirp_up=exp(1i*2*pi*fc_t.*time);       % chirp montant (utilise comme chirp brut)
chirp_down=exp(-1i*2*pi*fc_t.*time);    % chirp descendant


%% Emetteur
bits=rand(Ns,SF)>0.5;                   % genreation des bits
symb=bi2de(bits,'left-msb');            % generation des symboles
symb = gray2bin(symb,'pam',2^SF);       % passage en gray
Nb_preambule_up=5;                          % Nombre de chirp brut montant dans le preambule
Nb_preambule_down=2;                        % Nombre de chirp brut descendant dans le preambule

preambule=[repmat(chirp_up,1,Nb_preambule_up), repmat(chirp_down,1,Nb_preambule_down)]; % Signal d'apprentissage (header loRa)
sig=preambule;
PreambLegnth=length(preambule);

%% G?n?ration de l'enveloppe complexe ? partir de chirp montant
for k=1:length(symb)
    tau_k=symb(k)/B;
    time1=time(1:alpha*symb(k));
    fc_t1=B/(2*Ts)*(time1-2*tau_k)+B;
    time2=time(alpha*symb(k)+1:end);
    fc_t2=B/(2*Ts)*(time2-2*tau_k);
    sig=[sig exp(1i*2*pi*fc_t1.*time1) exp(1i*2*pi*fc_t2.*time2)];
end


