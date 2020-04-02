% Kelly Jawork, 4.2.20
% Need get3Dorbit, sv_from_oe, plot3Dorbit, plotSat, plotEverything
% To Do:
% - Add additional orbits & play with orbital parameters
% - Show cone of coverage, have side by side ground track
% - Satellite motion DONE
% - Shading on Mars ;)
% - Allow orbits to evolve over time due to perturbations
% - Fix legend
% - Keep axes fixed though figure is changing?
%% ORBITAL MANEUVERS USING LAMBERT'S PROBLEM 
clear;
clc;
close all;
%% USER INPUTS 
a = 8000; e = 0.02; % semimajor axis, eccentricity
i1 = 55; i2 = 110; i3 = 165; i4 = 215; % inclinations
RA = 60; w = 60; TA = 0; TA1 = 0; TA2 = 180; TA3 = 90; TA4 = 270;

for t=0:100 % iterates true anomaly over time
    % Plot axes, legend, Mars
    plotEverything()
    
    % Plot orbits
    plot3Dorbit(a,e,i1,RA,w,TA,'r-');
    plot3Dorbit(a,e,i2,RA,w,TA,'y-');
    plot3Dorbit(a,e,i3,RA,w,TA,'g-');
    plot3Dorbit(a,e,i4,RA,w,TA,'b-');
    
    % Plot satellites
    plotSat(a,e,i1,RA,w,TA1+t,'r-');
    plotSat(a,e,i2,RA,w,TA2+t,'y-');
    plotSat(a,e,i3,RA,w,TA3+t,'g-');
    plotSat(a,e,i4,RA,w,TA4+t,'b-');
    pause(0.000001);
    clf % clears figure; everything replotted with next time step
end
