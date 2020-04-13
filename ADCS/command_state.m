function u = command_state(e,e_sum,e_rate)
% y is a column vec:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16
% e: error
% e_sum: integrated error
% e_rate: error derivative

% 1) Gains
KP = 2; % proportional gain
KI = 0.01; % integral gain
KD = 0.01; % derivative gain

% 2) Command
u = KP*e + KI*e_sum + KD*e_rate; 
% u is a 3x1 vector
end

