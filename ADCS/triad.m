function q_triad = triad(N_sun, N_star, B_sun, B_star)
%{
Based on time and the location of the spacecraft, the sun and star tracker 
vectors in the reference frame N are provided. Based on the sun sensor and
star tracker readings, the sun and star tracker vectors are
provided in the spacecraft's body frame B. 

Because the TRIAD algorithm has the orientation of the same vectors in
different reference frames, it can compute the rotation between
those reference frames. The algorithm will fail if either of the reference
vectors or observed vectors are parallel or anti parallel.
%}

% Normalize all vectors (should be normalized already)
N_sun=N_sun/norm(N_sun);
B_sun=B_sun/norm(B_sun);
N_star=N_star/norm(N_star);
B_star=B_star/norm(B_star);

% Ensure all vectors have three components and vectors in the same frame
% are not parallel
if size(N_sun)~=[3 1]
    error('N_sun must be a 3-component column vector')
end

if size(B_sun)~=[3 1]
    error('B_sun must be a 3-component column vector')
end

if size(N_star)~=[3 1]
    error('N_mag must be a 3-component column vector')
end

if size(B_star)~=[3 1]
    error('B_mag must be a 3-component column vector')
end

if dot(N_sun,N_star)==1
    error('Reference frame N vectors cannot be parallel')
end

if dot(B_sun,B_star)==1
    error('Body frame B vectors cannot be parallel')
end

% Compute orthogonal right-handed vectors of reference frame
v1=N_sun;
v2=cross(N_sun,N_star)/norm(cross(N_sun,N_star));
v3=cross(v1,v2);

% Compute orthonogoal right-handed vectors of body frame
w1=B_sun;
w2=cross(B_sun,B_star)/norm(cross(B_sun,B_star));
w3=cross(w1,w2);

DCM_be=w1(:)*v1(:)'+w3(:)*v3(:)'+w2(:)*v2(:)';
q_triad=dcm2quat(DCM_be);
end
