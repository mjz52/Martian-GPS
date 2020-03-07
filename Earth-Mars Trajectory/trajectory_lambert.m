% Michael Zakoworotny
% MAE 4160
% Trajectory solver

function trajectory_lambert()

daySec = 3600*24; %seconds per day
% Sun properties
mu_S = 1.32712440018E20/(1000^3);   %G*Mass of sun (km^3/s^2)
R_S = 6.95700e5; %mean radius of Sun km

% Earth properties
mu_e = 398600.435436; %km^3/s^2
R_e = 6371.01; %mean radius of Earth in km
ae = 6378.137;          %km Earth geoid sma
fe = 1.0/298.257223563; %Earth geoid flattening
ee2 = 2*fe - fe^2;      % Earth geoid eccentricity squared
we = 7.292115e-5;       %rad/s rotation rate of the Earth

% Mars properties
mu_m = 42828.375214; %km^3/s^2
R_m = 3389.92; %km

global posE
global timeE
t_orbit_E = [datetime(2029,01,01,0,0,0), datetime(2031,12,31,0,0,0)]; %y,m,d,h,mi,s
% [timeE, posE, velE] = auto_horizons('EARTH', 'BARYCENTER', t_orbit_E(1), t_orbit_E(2));
% dlmwrite('posE.txt',posE); dlmwrite('velE.txt', velE); dlmwrite('timeE.txt', timeE);
posE = importdata('posE.txt'); velE = importdata('velE.txt'); timeE = importdata('timeE.txt');

global posM
t_orbit_M = t_orbit_E;
% [timeM, posM, velM] = auto_horizons('MARS', 'BARYCENTER', t_orbit_M(1), t_orbit_M(2));
% dlmwrite('posM.txt',posM); dlmwrite('velM.txt', velM);
posM = importdata('posM.txt'); velM = importdata('velM.txt');

% [timeP, posP, velP] = auto_horizons('Phobos', 'BARYCENTER', t_orbit_E(1), t_orbit_E(2));
% dlmwrite('posP.txt',posP); dlmwrite('velP.txt', velP);
% posP = importdata('posP.txt'); velP = importdata('velP.txt');

% IDENTIFY MINIMUM ENERGY TRANSFER
% t0_earliest = juliandate(t_orbit_E(1));
% t0_latest = juliandate(t_orbit_E(2));
% st = linspace(t0_earliest, t0_latest, 1000); %Array of start dates
% dur = @(t0) getMinOrbit(t0, mu_S); %Function that outputs duration given a start date
% %Get an array of durations of min energy orbit for each start date
% mindurs = zeros(length(st),1);
% for i = 1:365%length(st)
%     try
%         mindurs(i) = dur(st(i)); 
%     catch
%         mindurs(i) = length(st);
%     end
% end
% %Get minumum orbit
% minDur = min(mindurs);
% minSt = st(mindurs == minDur);
% minSt_date = datetime(minSt, 'ConvertFrom', 'juliandate');

fig = figure(1); shg;
[x_sun, y_sun, z_sun] = sphere();
x_sun = x_sun*R_S; y_sun = y_sun*R_S; z_sun = z_sun*R_S;
times = linspace(t_orbit_E(1),t_orbit_E(2),length(posE));
% Animate motion of Earth and Mars
% for i = 1:length(posE)
%     % The whole path length of the orbits
%     plot3(posE(1:365,1),posE(1:365,2),posE(1:365,3), 'LineWidth', 1, 'color', getColor('LightBlue'));
%     hold on;
%     plot3(posM(1:687,1),posM(1:687,2),posM(1:687,3), 'LineWidth', 1, 'color', getColor('Crimson'));
%     plot3(posP(1:687,1),posP(1:687,2),posP(1:687,3), 'LineWidth', 1, 'color', getColor('Grey'));
% 
%     % The current position of the planets
%     plot3(posE(i,1),posE(i,2),posE(i,3), '.', 'MarkerSize',20, 'color', getColor('Blue'));
%     plot3(posM(i,1),posM(i,2),posM(i,3), '.', 'MarkerSize',20, 'color', getColor('Red'));
%     plot3(posP(i,1),posP(i,2),posP(i,3), '.', 'MarkerSize',10, 'color', getColor('Grey'));
%     
%     % The sun
%     surf(x_sun, y_sun, z_sun);
%     
%     text(0,0,datestr(times(i)));
%     
%     hold off
%     axis equal;
%     % USE THIS TO ZOOM IN ON PLANETS
% %     axis([getBounds([posE(i,1),posM(i,1),posP(i,1)]) ...
% %             getBounds([posE(i,2),posM(i,2),posP(i,2)]) ...
% %             getBounds([posE(i,3),posM(i,3),posP(i,3)])]);
% 
% %     ylabel('$y$','interpreter','latex');
% %     xlabel('$x$','Interpreter','latex');
% %     title('Orbit of Two Particles in Space');
%     %legend({'Mass 1', 'Mass 2'}, 'Location', 'east');
%     perframe = 5/length(posE); %Seconds per frame to get a 10 second animation
%     pause(perframe); 
%     
%     %Ends loop if figure is closed
%     disp(~ishghandle(fig));
%     if ~ishghandle(fig)
%         break
%     end
% end

% GET ORBITAL POSITIONS
t1 = datetime(2026,11,08,0,0,0); %y,m,d,h,mi,s
t2 = datetime(2027,09,08,0,0,0);
% t1 = minSt_date;
% t2 = minSt_date + days(minDur);
% [time1, pos1, vel1] = auto_horizons('EARTH', 'BARYCENTER', t1, 0);
% [time2, pos2, vel2] = auto_horizons('MARS', 'BARYCENTER', t2, 0); 
pos1 = [104272285.062859,104470951.604389,6668.24206011742];
vel1 = [-21.6237048612658,20.8863036792556,-0.000919766639821518];
pos2 = [-101291740.609022,-203638179.164498,-1762748.03834854];
vel2 = [22.5963552736314,-8.73163333299846,-0.737188174778960];
time1 = 2461352.50000000;
time2 = 2461656.50000000;

v_orbit_E = sqrt(9 + 2*mu_e/(R_e+200));
v1_ = vel1' - vel1'/norm(vel1')*v_orbit_E;
v_orbit_M = sqrt(0 + 2*mu_m/(R_m+20000));
v2_ = vel2' - vel2'/norm(vel2')*v_orbit_M;
[vi, vf] = glambert(mu_S, [pos1, vel1], [pos2, vel2], (time2-time1)*daySec, 0);
% [vi, vf,extremal_distances, exitflag] = lambert(pos1,pos2,(time2-time1),0,mu_S)

[~, z_array] = ode45(@newtongravity2d,[time1*daySec, time2*daySec],[pos1, vi'], ...
                    odeset('RelTol',1e-12,'AbsTol',1e-12), mu_S);
xTR1 = z_array(:,1); yTR1 = z_array(:,2); zTR1 = z_array(:,3);

figure(2); hold on; axis equal;
plot3(posE(:,1),posE(:,2),posE(:,3),'color','blue');
plot3(posM(:,1),posM(:,2),posM(:,3),'color','red');
plot3(pos1(1),pos1(2),pos1(3),'*b');
plot3(pos2(1),pos2(2),pos2(3),'*r');
plot3(xTR1,yTR1,zTR1,'color', 'black');



delv1 = norm(vi - v1_);
delv2 = norm(v2_ - vf);
fprintf('Delta v injection is %.2f km/s\n',delv1);
fprintf('Delta v arrival is %.2f km/s\n',delv2);



%% Orbital propagator
    function dz = newtongravity2d(~, z, mu)
        %Input: state vector z
        %output: state derivative dz
        d = sqrt(z(1)^2 + z(2)^2 + z(3)^2);
        ax = -mu/d^3 * z(1);
        ay = -mu/d^3 * z(2);
        az = -mu/d^3 * z(3);
        dz = [z(4), z(5), z(6), ax, ay, az].';
    end

%% Minimum energy orbit helper functions
    function dur = getMinOrbit(st, mu)
        %Given an array of starting times, finds the "zeros" of the function,
        %ie. trajectory durations that satisfy the equations
        minTime = @(dt) getOrbit(st,dt,mu) - dt;
        dur = fzero(minTime, 0); %Search around 
    end


    function tm = getOrbit(st, dur, mu)
        %Output the minimum transfer time given conditions
        [r1_, r2_] = EarthMarsPositions(st, dur); %Get r1 and r2
    %         orb.r1_ = r1_*p.kmAU; orb.r2_ = r2_*p.kmAU;
        r1 = norm(r1_); r2 = norm(r2_); c = norm(r1_ - r2_);
        s = (r1+r2+c)/2;
        taum = pi*sqrt(s^3/(2*mu));
        betam = 2*asin(sqrt((s-c)/s));
        tm = taum/(2*pi)*(pi - (betam - sin(betam)));
    end


    function [r1_, r2_] = EarthMarsPositions(st, dur)
        ind = find(timeE == st);

        r1_ = posE(ind,:); r2_ = posM(ind+max(0,round(dur)),:);
    end


    function v = get_vel_circ(r,mu)
       v = sqrt(mu*2/r); 
    end


    function minmax = getBounds(vals)
        max = vals(1); min = vals(1); 
        for k = 1:length(vals)
            if (vals(k)<min)
               min = vals(k); 
            end
            if (vals(k)>max)
               max = vals(k); 
            end
        end
        diff = abs(min - max);
        minmax = [min-0.2*diff max+0.2*diff];
%         minmax = [sign(min)*0.8*abs(min), sign(max)*1.2*abs(max)];
    end

end