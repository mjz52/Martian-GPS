% Michael Zakoworotny
% MAE 4160

close all; clear; %clc;

%%
daySec = 3600*24; %seconds per day
J2000 = datetime(2000,01,01,12,0,0);

% Sun properties
mu_S = 1.32712440018e20/(1000^3);   %G*Mass of sun (km^3/s^2)
R_S = 6.95700e5; %mean radius of Sun km

% Earth properties
mu_e = 398600.435436; %km^3/s^2
R_e = 6371.01; %mean radius of Earth in km
ae = 6378.137;          %km Earth geoid sma
fe = 1.0/298.257223563; %Earth geoid flattening
ee2 = 2*fe - fe^2;      % Earth geoid eccentricity squared
we = 7.292115e-5;       %rad/s rotation rate of the Earth
be = ae*sqrt(1-ee2);
J2_e = 1082.645e-6;


% Mars properties
mu_m = 42828.375214; %km^3/s^2
R_m = 3389.92; %km
fm = 1/169.779; %Flattening
am = 3396.19; %Equitorial radius (km)
em2 = 2*fm - fm^2; %Geoid eccentricity squared
wm = 0.0000708822; %rotation rate (rad/s)
bm = am*sqrt(1-em2); %Polar radius (km)
C20 = -0.8750220924537000E-03; l = 2; m = 0;
J2_m = -sqrt(factorial(l-m)*(2*l+1)*(2-1)/factorial(l+m))*C20;

figure(1); hold on; axis equal;
%[pos,vel] = planetEphemeris(EPHEMERISTIME, CENTER, TARGET, MODEL, UNITS, ACTION)
% CORRESPONDING SETTINGS IN HORIZONS:
% Reference system: ICRF/J2000
% Reference X-Y plane: Earth Mean Equator and Equinox of Reference Epoch
% 405 - J200 from 1997
% 430 - most up-to-date ephemerides
% Plot Earth orbit
% t1 = datetime(2020,01,01,00,00,00);
% t2 = t1+years(5);
% t1E = datetime(2017,09,12,0,0,0); %EXAMPLE FROM ONLINE PICTURE
% t2E = datetime(2018,10,17,0,0,0);
t1E = datetime(2028,01,01,0,0,0);
t2E = datetime(2036,01,01,0,0,0);
t_e = [t1E:t2E]';%linspace(t1,t2,366)';
[posE,velE] = planetEphemeris(juliandate(t_e),'Sun','Earth','430','km');
% posE = -posE; velE = -velE; %NEED TO NEGATE VECTOR
plot3(posE(:,1),posE(:,2),posE(:,3));

cent_E = (juliandate(mean(t_e))-juliandate(J2000))/(365.24*100);
posE = equitorial2ecliptic(posE,cent_E);
velE = equitorial2ecliptic(velE,cent_E);
x_e = posE(:,1); y_e = posE(:,2); z_e = posE(:,3);
plot3(x_e,y_e,z_e);  % Plot orbits in ecliptic
i = find(t_e == datetime (2020,03,19,0,0,0));
text(posE(i,1),posE(i,2),posE(i,3),"Mar-19",'FontSize',8);
i = find(t_e == datetime (2020,06,20,0,0,0));
text(posE(i,1),posE(i,2),posE(i,3),"Jun-20",'FontSize',8);
i = find(t_e == datetime (2020,09,22,0,0,0));
text(posE(i,1),posE(i,2),posE(i,3),"Sep-22",'FontSize',8);
i = find(t_e == datetime (2020,12,21,0,0,0));
text(posE(i,1),posE(i,2),posE(i,3),"Dec-21",'FontSize',8);

% Plot Mars orbit
% t1 = datetime(2020,01,01,00,00,00);
% t2 = t1+years(7);
% t1M = datetime(2018,4,15,0,0,0); %EXAMPLE FROM ONLINE PICTURE
% t2M = datetime(2020,3,15,0,0,0);
t1M = datetime(2028,06,01,0,0,0);
t2M = datetime(2037,06,01,0,0,0);
t_m = [t1M:t2M]';%linspace(t1,t2,688)';
[posM,velM] = planetEphemeris(juliandate(t_m),'Sun','Mars','430','km');
% posM = -posM; velM = -velM; %NEED TO NEGATE VECTOR

cent_M = (juliandate(mean(t_m))-juliandate(J2000))/(365.24*100);
posM = equitorial2ecliptic(posM,cent_M);
velM = equitorial2ecliptic(velM,cent_M);
x_m = posM(:,1); y_m = posM(:,2); z_m = posM(:,3);
plot3(x_m,y_m,z_m);

drawAxes();
i = find(y_e==max(y_e));
max_pos = [x_e(i),y_e(i),z_e(i)];
inc = acosd((max_pos(2)/norm(max_pos)));
j = find(x_e==max(x_e));
equinox = [x_e(j),y_e(j),z_e(j)];

dv_max = 15.5;

pork = zeros(length(t_e),length(t_m));
del_t = days(t1M-t1E);
for i = 1:length(t_e)
    i
    t_aft = find(t_m>t_e(i),1);
    for j = t_aft:length(t_m)
        sv1 = [posE(i,:), velE(i,:)];
        sv2 = [posM(j,:), velM(j,:)];
        try
            [vi,vf] = glambert(mu_S, sv1, sv2, (j-i+del_t)*daySec, 0);
            dv1 = norm(sv1(4:6) - vi');
            dv2 = norm(vf' - sv2(4:6));
            dv = dv1 + dv2;
            pork(i,j) = (dv<dv_max)*dv;
        catch
            pork(i,j) = 0;
        end
    end
end

%% Plotting porkchop
[T_m, T_e] = meshgrid(t_m,t_e);
figure;
T_e = days(t_e - datetime(0,1,0,0,0,0));
T_m = days(t_m - datetime(0,1,0,0,0,0));
x = repelem((T_e(1):T_e(end))',1,length(t_m));
y = repelem((T_m(1):T_m(end)),length(t_e),1);


[M,c] = contourf(x,y,pork);

% co = [getColor('LightBlue');getColor('Aqua');getColor('Navy');getColor('DarkBlue')];
set(gca,'XLim',[T_e(1) T_e(end)]);
set(gca,'YLim',[T_m(1) T_m(end)]);
set(gca,'XTick',ceil(linspace(T_e(1),T_e(end),round(years(t_e(end)-t_e(1)))+1)));
set(gca,'YTick',(linspace(T_m(1),T_m(end),round(years(t_m(end)-t_m(1)))+1)));
datetick('x','dd-mmm-yy','keeplimits','keepticks'); datetick('y','dd-mmm-yy','keeplimits','keepticks');
set(gca,'TickDir','out');

xlabel('Departure Date (from Earth)');
ylabel('Arrival Date (at Mars)');
title('Porkchop Plot for Earth and Mars with \Deltav contours', 'interpreter', 'tex', 'fontname', 'Ti89pc')
colorbar();

%% Find local minimima
val = []; r = []; c = [];
tol = days(20);
% Get regions of time when there are valid departure and arrival dates
delv_acceptable = 7.5;%round(min(min(pork(pork~=0))),1); %set minimum to lowest attainable del_v
while ~(length(val) > 10 && length(find(diff(sort(t_e(r)))>tol)) >= 6)
    delv_acceptable = delv_acceptable + 0.1
%     TF = islocalmin(pork,1) & islocalmin(pork,2) & pork < delv_acceptable & pork~=0;
    TF = islocalmin(pork,1) & islocalmin(pork,2) & pork < delv_acceptable & pork~=0;
    [r,c,val] = find(TF);
end

[t_e_sort, i_sort] = sort(t_e(r));
sections = [0; find(diff(t_e_sort)>tol); length(t_e_sort)];
dep_dates = t_e(r); dep_dates = dep_dates(i_sort);
arr_dates = t_m(c); arr_dates = arr_dates(i_sort);
del_v = zeros(length(r),1);
for i = 1:length(r)
	del_v(i) = pork(r(i),c(i));
end
del_v = del_v(i_sort);
% Get best times from each section
best_dep_dates = [];
best_arr_dates = [];
best_del_v = [];
for j = 1:length(sections)-1
    del_v_sect = del_v(sections(j)+1:sections(j+1));
    [~, ind_sort] = sortrows(del_v_sect);
    k = ind_sort(1:min(3,length(ind_sort)));
    k = k + sections(j);
    best_dep_dates = [best_dep_dates; dep_dates(k)];
    best_arr_dates = [best_arr_dates; arr_dates(k)];
    best_del_v = [best_del_v; del_v(k)];
end

best_dates = [best_dep_dates, best_arr_dates]
best_del_v

%% Plot each date
hold on;
for i = 1:length(best_dates)
    td = find(t_e==best_dep_dates(i));
    ta = find(t_m==best_arr_dates(i));
    plot(T_e(td),T_m(ta),'*','color','black')
end



%% Helper functions

function drawAxes()
    k = 1.1;
    line([0 k*max(xlim)], [0 0], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 k*max(ylim)], [0 0], 'LineWidth', 3, 'color', 'black');
    line([0 0], [0 0], [0 k*max(zlim)], 'LineWidth', 3, 'color', 'black');
    text(k*max(xlim),0,0,"X-axis",'FontSize',8);
    text(0,k*max(ylim),0,"Y-axis",'FontSize',8);
    text(0,0,k*max(zlim),"Z-axis",'FontSize',8);
%     leg = legend(gca);
%     leg = leg.String;
%     leg(end-2:end) = {"X-axis","Y-axis","Z-axis"};
%     legend(leg, 'Location', 'eastoutside');
end

