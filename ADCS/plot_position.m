function plot_position(C,t,state_arr,title_str)
% state_arr has:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16
x = state_arr(1,:)'; y = state_arr(2,:)'; z = state_arr(3,:)';
figure()
subplot(3,1,1)
plot(t,x,'b--')
title(title_str)
ylabel('x')
subplot(3,1,2)
plot(t,y,'r--')
ylabel('y')
subplot(3,1,3)
plot(t,z,'g--')
ylabel('z')
end
