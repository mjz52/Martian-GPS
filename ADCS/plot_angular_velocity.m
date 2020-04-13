function plot_angular_velocity(C,t,state_arr)
% state_arr has:
% [x y z vx vy vz q1 q2 q3 q4 wx wy wz wgx wgy wgz]
%  1 2 3 4  5  6  7  8  9  10 11 12 13 14  15  16
wx = state_arr(11,:)'; wy = state_arr(12,:)'; wz = state_arr(13,:)';
figure()
subplot(3,1,1)
plot(t,wx,'b--')
ylabel('\omega_1')
subplot(3,1,2)
plot(t,wy,'r--')
ylabel('\omega_2')
subplot(3,1,3)
plot(t,wz,'g--')
ylabel('\omega_3')
end

