% Rotate a vector by the inclination. Uses definition of obliquity of
% inclination from Astronomical Almanac 2010
% r: N x 3 matrix containing position vector coordinates (alternatively
% velocity)
% T: time period, in centuries since J2000
function r_ecl = equitorial2ecliptic(r,T)
    % Angle of inclination of 
    eps = dms2degrees([23 26 21.406]) - dms2degrees([0 0 46.836769])*T - ...
          dms2degrees([0 0 0.0001831])*T.^2 + dms2degrees([0 0 0.00200340])*T.^3 - ...
          dms2degrees([0 0 0.576e-6])*T.^4 - dms2degrees([0 0 4.34e-8])*T.^5;
    C = [1, 0, 0;
         0, cosd(eps), sind(eps);
         0, -sind(eps), cosd(eps)];
    r_ecl = C*r';
    r_ecl = r_ecl';
end