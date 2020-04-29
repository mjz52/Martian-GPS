
close all; clear; clc;
%% Power Calculations

Eb_N0 = 5; %Recommended
Gt = 10.3+10*log10(52/10); %dB
Gr = 5; % 3 to 5 dB
La = -1; %dB
Ll = -2; %dB
Lamb = 3e8/1575.42e6; %MHz
Kb = 1.38064852e-23; % m^2 kg /s^2 K
Tsys = 243; %K, Marian Atmosphere ranges from 120K to 200K
Rb = 1.02389e6; %bps
d = linspace(0,50000e3,100)'; del_d = diff(d); del_d = del_d(1);
Pt = Eb_N0 - Gt - Gr - La - Ll - 10*log10((Lamb./(4*pi*d)).^2)+10*log10(Kb*Tsys*Rb);
Pt = 10.^(Pt/10);
d_power = d;

%% Satellite Quantity Values
satQty = [11 12 13 14 15 16 17];
d_sat = 1e3*[50000 22000 13000 11000 7000 6000 5000];

%% Intersection of two plots
satQty_scale = (satQty-min(satQty)) *(max(Pt)-min(Pt))/(max(satQty)-min(satQty));
x = [flip(d_sat(2:3))', flip(satQty_scale(2:3))'];
y = [d_power(43:44), Pt(43:44)]
% y = [21210000, 2424;
%      21720000, 2541];
X = fzero(@(X)(y(2,2)-y(1,2))/(y(2,1)-y(1,1))*(X-y(1,1))+y(1,2)-((x(2,2)-x(1,2))/(x(2,1)-x(1,1))*(X-x(1,1))+x(1,2)), 2e7);
di = X;
pi_ = interp1(y(:,1),y(:,2),X);

%% Plotting
% Plot altitude vs # satellites along with 
fig = figure; hold on;

yyaxis left
ax1 = plot(d_sat,satQty,'-o','LineWidth',2,'color','b');
ylabel('Number of satellites');
yyaxis right
ax2 = plot(d,Pt,'--','LineWidth',2,'color','r');
plot(di,pi_,'^','MarkerSize',8,'color','black');
text(2.4e7,250,['N_{sat} = ' num2str(12) ', P = ' num2str(pi_,'%.1f') ' W'],'interpreter','tex')
ylabel('Power [W]');
xlabel('Altitude [m]');
% title('Satellite quantity and power vs altitude');
legend({'Number of satellites','Power','Selected Configuration'}, 'location','N');
ax = gca;
set(ax.YAxis(1),'color','black')
set(ax.YAxis(2),'color','black')
set(ax.YAxis(1),'Limits',[min(satQty),max(satQty)]);
set(ax.YAxis(2),'Limits',[min(Pt),max(Pt)]);
set(ax,'FontSize',11)
set(ax,'FontName','Ti89pc');

%% Plot # Satellites vs Power
Pt_sat = zeros(1,length(satQty))
ind = []
for i = 1:length(d_sat)
    ind = [ind find(abs(d-d_sat(i))<del_d/2)];
end
Pt_sat = Pt(ind)
figure;
plot(satQty,Pt_sat);
% title('Satellites power vs Satellites quantity');

%% Plot product
% plot(d_sat, Pt_sat.*satQty')


%% Antenna Sizing
z = linspace(0.8,1.2,100);
% L = (27.04*Lamb)./(z.^2);
D = (Lamb.*z)/pi;
L = Lamb^3*(52./D/pi/10).^2;
figure; hold on;
plot(D,L,'LineWidth',2,'color','black');
plot(max(D),min(L),'^','MarkerSize',8,'color','black');
text(0.067,4.4,['D = ' num2str(max(D),'%.3f') 'm, L = ' num2str(min(L),'%.1f') ' m'],'interpreter','tex')
xlabel('D [m]');
ylabel('L [m]');
% title('Helical Antenna Sizing for G_T=17.46 and \alpha=10^{\circ}','interpreter','tex');
legend({'Possible configurations','Selected configuration'},'location','NE');
ax = gca;
% set(ax.YAxis(1),'Limits',[min(L),max(L)]);
% set(ax.YAxis(2),'Limits',[min(D),max(D)]);
set(ax,'FontSize',11)
set(ax,'FontName','Ti89pc');

% plot(z,L)
% plot(z,D)
% hold off
%% Intersection of two plots
L_scale = (L-min(L))*(max(D)-min(D))/(max(L)-min(L)) + min(D)
x = [flip(z(43:44))', flip(L_scale(43:44))'];
y = [flip(z(43:44))', flip(D(43:44))'];
% y = [0.9697, 7.591e-5;
%      0.9737, 7.623e-5];
X = fzero(@(X)(y(2,2)-y(1,2))/(y(2,1)-y(1,1))*(X-y(1,1))+y(1,2)-((x(2,2)-x(1,2))/(x(2,1)-x(1,1))*(X-x(1,1))+x(1,2)), 2e7);
zi = X;
Di = interp1(y(:,1),y(:,2),X)
Li = interp1(y(:,1),L(43:44)',X)

%% Antenna Plot
figure; hold on;
yyaxis left
ax1 = plot(z,L,'-.','LineWidth',2,'color','b');
ylabel('L [m]');
yyaxis right
ax2 = plot(z,D,'--','LineWidth',2,'color','r');
plot(zi,Di,'^','MarkerSize',8,'color','black');
text(1,0.059,['L = ' num2str(Li,2) ' m, D = ' num2str(Di,2) ' m'])
ylabel('D [m]');
xlabel('z');
% title('Antenna Sizing');
legend({'L','D','Intersection'}, 'location','N');
ax = gca;
set(ax.YAxis(1),'color','black')
set(ax.YAxis(2),'color','black')
set(ax.YAxis(1),'Limits',[min(L),max(L)]);
set(ax.YAxis(2),'Limits',[min(D),max(D)]);
set(ax,'FontSize',11)
set(ax,'FontName','Ti89pc');


% figure
% Hif = linspace(60,100,1000);
% G = Gt + 20;
% HiD = 10.^((1/20)*(G-17.8-20*log10(Hif)));
% plot(Hif,HiD)
%%
Hif = linspace(60,100,1000);
HiG = Gt+20;
HiD = 10.^((1/20)*(HiG-17.8-20*log10(Hif)));
figure;
plot(Hif,HiD,'--','Linewidth',2,'color','black');
xlabel('\lambda [GHz]','interpreter','tex');
ylabel('D [m]','interpreter','tex');
ax = gca;
set(ax,'FontSize',11)
set(ax,'FontName','Ti89pc');


Eb_N0 = 5; %Recommended
Gr = 5; % 3 to 5 dB
La = -1; %dB
Ll = -2; %dB
Lamb = 3e8/1575.42e6; %m
Kb = 1.38064852e-23; % m^2 kg /s^2 K
Tsys = 200; %K, Marian Atmosphere ranges from 120K to 200K
Rb = 448*10000; %bps
d = 2.5e7;
Pt = Eb_N0 - Gt - Gr - La - Ll - 10*log10((Lamb/(4*pi*d))^2)+10*log10(Kb*Tsys*Rb);
Pt = 10.^(Pt/10);

