function [tarray,zarray] = ode_solver(f,state0,sensors,actuators)
% Runge-Kutta 4th order method
global const

dt = const.dt;
h = dt;  % set the step size
x = 0:dt:const.TEND;  % set the interval of x
y = zeros(1,length(x));
y(0) = state0;   % set the intial value for y
n = length(x)-1;
y_dot =@(x,y)(f); %insert function to be solved
for i = 1:n
    
    k1 = y_dot(x(i),y(i));
    k2 = y_dot(x(i)+.5*h,y(i)+.5*k1*h);
    k3 = y_dot(x(i)+.5*h,y(i)+.5*k2*h);
    k4 = y_dot(x(i)+h,y(i)+k3*h);
    y(i+1) = y(i)+((k1+2*k2+2*k3+k4)/6)*h;
    
end
end

