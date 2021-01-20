function d = calcul_de_la_distance(f,S,P_tx)

c = 3*10^8;

lambda = c/f;

PL = - P_tx - S;

d = lambda/(4*pi*10^(PL/20));

end