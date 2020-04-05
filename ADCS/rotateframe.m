function u = rotateframe(q,v)
% u = rotateframe(q,v): rotates the frame of reference for the (3x1) 
% vector v using quaternion q stored as an array with the 4th component real.

u = v  + cross(2*q(1:3), cross(q(1:3),v) - q(4)*v);
end