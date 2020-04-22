Eb_N0 = 5; %Recommended
Gt = 13; %dB
Gr = 5; % 3 to 5 dB
La = -1; %dB
Ll = -2; %dB
Lamb = 3e8/1575.42e6; %MHz
Kb = 1.38064852e-23; % m^2 kg /s^2 K
Tsys = 200; %K, Marian Atmosphere ranges from 120K to 200K
Rb = 448*10000; %bps

d = linspace(0,50000e3,1000)'; del_d = diff(d); del_d = del_d(1);
Pt = Eb_N0 - Gt - Gr - La - Ll - 10*log10((Lamb./(4*pi*d)).^2)+10*log10(Kb*Tsys*Rb);
Pt = 10.^(Pt/10);
%plot(d,Pt)

% Plot altitude vs # satellites along with 
figure; hold on;
SatQty = [11 12 13 14 15 16 17];
min_distance = 1e3*[50000 22000 13000 11000 7000 6000 5000];
yyaxis left
plot(min_distance,SatQty);
ylabel('Number of satellites');
yyaxis right
plot(d,Pt);
ylabel('Power [W]')
title('Satellite quantity and power vs altitude');
legend({'Quantity of satellites','Satellite power'});

% Plot # Satellites vs Power
Pt_sat = zeros(1,length(SatQty))
ind = []
for i = 1:length(min_distance)
    ind = [ind find(abs(d-min_distance(i))<del_d/2)];
end
Pt_sat = Pt(ind)
figure;
plot(SatQty,Pt_sat);
title('Satellites power vs Satellites quantity');


