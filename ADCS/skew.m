function sk = skew(x)
% returns the skew symmetric matrix of x: [x1;x2;x3]
sk = [0 , -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
end

