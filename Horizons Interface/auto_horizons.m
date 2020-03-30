% Michael Zakoworotny
% 
% 

% body - String specifying the target body
% observer - String specifying the observer body
% st_time - starting time
% end_time - ending time. If 0, default to using a time minute interval
% (only want one state vector, not the entire list)

function [time, pos, vel] = auto_horizons(body, observer, st_time, end_time)

% st_time = juliandate(datetime(datestr(now)) + days(0)); %Get current time
st_time = juliandate(st_time);
try (end_time == 0);
    dur = minutes(10);
    end_time = juliandate(datetime (st_time,'ConvertFrom','juliandate') + dur);
    ret_all = 0; %Only return the first state
catch
    end_time = juliandate(end_time);
    ret_all = 1; %Return all states
end

file_name = 'auto_horizons.txt';

gen_file(body, observer, st_time, end_time, file_name, ret_all);

t=tcpip('horizons.jpl.nasa.gov',6775);
buffSize = 1000000;
set(t,'InputBufferSize',buffSize);
set(t,'OutputBufferSize',buffSize);
fopen(t);
pause(0.1);
fid = fopen(file_name);
while(~feof(fid))
    if(t.BytesAvailable > 0)
        text = fread(t,[1,t.BytesAvailable]);
        native2unicode(text);
        line = fgetl(fid);
        fprintf(t, line);
%         disp('==============================');
    end
end
if(ret_all == 0)
    pause(0.5);
else
    pause(10);
end
output = native2unicode(fread(t,[1,t.BytesAvailable]));
results = get_state(output);
time = results(:,1);
position = results(:,2:4);
velocity = results(:,5:7);

if (ret_all == 0)
    time = time(1,:); pos = position(1,:); vel = velocity(1,:);
else
    time = time(:,:); pos = position(:,:); vel = velocity(:,:);
end

end


%% GENERATE FILE TO RUN HORIZONS
function gen_file(target, observer, st_time, end_time, file_name, ret_all)
    
    fid = fopen(file_name,'w');
    fprintf(fid, 'major-bodies\n\n');
    
    t_id = get_target(target); % Target body
    fprintf(fid, [t_id '\n']);
    
    comp = 'v'; % vectors (set to e for orbital elements)
    fprintf(fid, ['E\n' comp '\n']); %Ephemerides data, vector components
    
    o_id = get_observer(observer); % Observer body
    fprintf(fid, [o_id '\ny\n']);
    
    plane = 'eclip'; % Or set to frame or body
    fprintf(fid, [plane '\n']);
    
    st = datestr(datetime(st_time, 'ConvertFrom', 'juliandate'));
    fprintf(fid, [st '\n']);
%     slice = regexp(st, ':');
%     st = st(1:slice(end)-1);
%     min = st(end-1:end);
%     hour = st(end-4:end-3);
%     new_min = min+dur;
%     if (new_min)
    et = datestr(datetime(end_time, 'ConvertFrom', 'juliandate'));
    fprintf(fid, [et '\n']);
    
    if (ret_all == 0)
        time_int = '10m';
    else
        time_int = '1d';
    end
    
    fprintf(fid, [time_int '\n']); % Polling frequency, set to '1h', '1d', etc
    fprintf(fid, 'n\nJ2000\n1\n'); %Reference frame, corrections
    
    units = '1'; % 1=km/s, 2=AU/D, 3=KM/D
    fprintf(fid, [units '\n']);
    fprintf(fid, 'YES\nYES\n'); % Output CSV format, Output delta-T
    
    output_code = '2'; %1=position, 2=state vector, 3=state vector+extra,
                       % 4=position+extra, 5=velocity only, 6=extra
                       % extra = downleg light time, range, range-rate
    fprintf(fid, [output_code '\n']);
%     fprintf(fid, '\nM\nmjz52@cornell.edu\nyes\n');
%     fprintf(fid, '\nq');
    
    fclose(fid);       
end

%% GET STATE VECTOR
% Parse state vector from Horizons info
function results = get_state(output)
    
    bounds = find(output == '$');
    data = output(bounds(2):bounds(3));
    returns = regexp(data, '[\r]');
    lines = {};
    for i = 1:length(returns)-1
        lines{i} = data(returns(i)+2:returns(i+1)-1);
    end
    
    results = zeros(length(lines),7);
    
    for i = 1:length(lines)
        %Get time
        nums = regexp(lines{i}, '(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?', 'match');
        time = str2num(nums{1});
        results(i,1) = time;
        %Get positions
        x = str2num(nums{8});
        y = str2num(nums{9});
        z = str2num(nums{10});
        results(i,2:4) = [x,y,z];
        %Get velocities
        vx = str2num(nums{11});
        vy = str2num(nums{12});
        vz = str2num(nums{13});
        results(i,5:7) = [vx,vy,vz];
    end
end


%% GET ID's
% Get ID of target body (ie. body for which you want Ephemeride)
function id = get_target(target)
    id = '0'; %Default to solar system Barycenter
    switch target
        case 'MERCURY'
            id = '199';
        case 'VENUS'
            id = '299';
        case 'EARTH'
            id = '399';
        case 'MARS'
            id = '499';
        case 'JUPITER'
            id = '599';
        case 'SATURN'
            id = '699';
        case 'URANUS'
            id = '799';
        case 'NEPTUNE'
            id = '899';
        case 'SUN'
            id = '10';
        case 'MOON'
            id = '301';
        case 'Phobos'
            id = '401';
        case 'Deimos'
            id = '402';
    end
end

% Get ID of observer location
function id = get_observer(obs)
    id = '@0';
    switch obs
        case 'BARYCENTER' % Solar system barycenter
            id = '@0';
        case 'SUN' % Center of sun
            id = '@sun';
        case 'EARTH'
            id = '500' % Earth geocentric
        case 'MARS'
            id = '500@4' % Mars barycenter
    end

end