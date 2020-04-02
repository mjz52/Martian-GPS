% Function that plots 3D orbit and satellite from input orbital parameters
% Kelly Jawork, 3.30.20
function plot3Dorbit(a,e,incl,RA,w,TA,color)
hold on
[r0,rx0,ry0,rz0] = get3Dorbit(a,e,incl,RA,w,TA);
plot3(rx0,ry0,rz0,color,'LineWidth', 1);
end
