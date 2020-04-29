% Battery sizing

DOD = 0.5;
n = 0.9;

mu =  42828.375214;
a = 28390;
T = 2*pi*sqrt(a^3/mu);
Te = T/2;

N_cycles = 30*(365.25)/(T/3600/24)

Pe = 31.5;

Cr = Pe*Te/DOD/n;
Cr_Whr = Cr/3600
Spec_Energy = 60; % Nickel-hygrogen
M_bat = Cr_Whr/Spec_Energy
