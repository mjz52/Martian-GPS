function plotSat(a,e,incl,RA,w,TA,color)
    [r,rx,ry,rz] = get3Dorbit(a,e,incl,RA,w,TA);
    % plot motion of satellite positions
    plot3(r(1),r(2),r(3),'ob', 'MarkerSize',8,'MarkerFaceColor','c');
end