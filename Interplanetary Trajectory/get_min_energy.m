% Michael Zakoworotny
% Get minimumum energy interplanetary orbit

mu = 1.32712440018e20/(1000^3);   %G*Mass of sun (km^3/s^2)

t0_earliest = juliandate(datetime(2020,01,01,00,00,00));
time_span = 365*20;
t0_latest = t0_earliest + time_span;
st = linspace(t0_earliest, t0_latest, time_span); %Array of start dates
dur = @(t0) getMinOrbit(t0, mu); %Function that outputs duration given a start date
%Get an array of durations of min energy orbit for each start date
mindurs = zeros(length(st),1);
for i = 1:length(st)
   mindurs(i) = dur(st(i));
end
%Get minumum orbit
minDur = min(mindurs);
minSt = st(mindurs == minDur);
minSt_date = datetime(minSt, 'ConvertFrom', 'juliandate');



function dur = getMinOrbit(st, mu)
    %Given an array of starting times, finds the "zeros" of the function,
    %ie. trajectory durations that satisfy the equations
    minTime = @(dt) getOrbit(st,dt,mu) - dt;
    dur = fzero(minTime, 0); %Search around 
end

function tm = getOrbit(st, dur, mu)
    %Output the minimum transfer time given conditions
    st
    dur
    r1_ = planetEphemeris(st,'Sun','Earth','430','km');
    r2_ = planetEphemeris(st+dur,'Sun','Mars','430','km');
    r1 = norm(r1_); r2 = norm(r2_); c = norm(r1_ - r2_);
    s = (r1+r2+c)/2;
    taum = pi*sqrt(s^3/(2*mu));
    betam = 2*asin(sqrt((s-c)/s));
    tm = taum/(2*pi)*(pi - (betam - sin(betam)));
end