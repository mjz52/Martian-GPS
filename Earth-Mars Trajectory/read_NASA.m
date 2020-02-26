% Return a matrix containing the position and veloctiy vectors stored in a
% JPL Horizons data file.
% TODO:
% - Implement the type parameter to differentiate between files containing
% orbital elements and orbital vectors. Duplicate the current orbital
% vector code for the orbital elements case

function results = read_NASA(file, type)

% addpath('C:\Users\mjz17\Documents\Michael\Semester 5\Spaceflight\Project');
fid = fopen(file);
k = 1; data = {};
bounds = [];
while(~feof(fid))
    data{k} = fgetl(fid);
    if(length(data{k})>0 && data{k}(1) == '$')
        bounds = [bounds; k];
    end
    k = k+1;
end
data = data';

results = zeros(mod(bounds(2)-bounds(1),4),7);
k = 1;
for i = bounds(1)+1 : 4 : bounds(2)-1
    %Get time
    nums = regexp(data{i}, '(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?', 'match');
    time = str2num(nums{1});
    results(k,1) = time;
    %Get positions
    nums = regexp(data{i+1}, '(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?', 'match');
    x = str2num(nums{1});
    y = str2num(nums{2});
    z = str2num(nums{3});
    results(k,2:4) = [x,y,z];
    %Get velocities
    nums = regexp(data{i+2}, '(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?', 'match');
    vx = str2num(nums{1});
    vy = str2num(nums{2});
    vz = str2num(nums{3});
    results(k,5:7) = [vx,vy,vz];
    k = k+1;
end


