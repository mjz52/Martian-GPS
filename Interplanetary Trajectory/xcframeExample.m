
function xcframeExample(pos, time)
% xcframeExample(pos, time)
%
% This is an example on how to use the xcframe function to convert position
% between different coordinate frames. Optinonal input "pos" is a position on
% earth given as a three-element vector with decimal latitude/longitude in
% degrees and altitude in meters, default is the position for London, UK.
% Optional input time is a 6 element time vector [year month day hour minutes
% seconds] with UTC time, default is the current time.
% Default position
if (nargin < 1 || isempty(pos))
  % London, UK
  pos =  [51.509865, -0.118092, 0];
end
% Default time
if (nargin < 2 || isempty(time))
  % Current UTC time
  time = datevec(datetime('now', 'TimeZone', 'UTC'));
end
% Convert input position to spherical and Cartesian coordinates
sphPos  = xcframe([], pos, [], [], 'GEO', 'GEO', 'geod84', 'sph' );
cartPos = xcframe([], pos, [], [], 'GEO', 'GEO', 'geod84', 'cart');
% Display input position
clc
disp(['Position in the GEO frame at ', datestr(time)])
% disp(['Geodetic coordinates: [', num2str(pos), ']'])
fprintf('Geodetic coordinates : [%.2f, %.2f, %.2f] (deg, deg, km)\n', ...
  pos(1), pos(2), pos(3)*1e-3)
fprintf('Spherical coordinates: [%.2f, %.2f, %.2f] (deg, deg, km)\n', ...
  sphPos(1), sphPos(2), sphPos(3)*1e-3)
fprintf('Cartesian coordinates: [%.2f, %.2f, %.2f] (km)\n', ...
  cartPos(1)*1e-3, cartPos(2)*1e-3, cartPos(3)*1e-3)
% A cell array with all supported coordinate frames to convert to
frames = {'GEI', 'J2000', 'MAG', 'GSE', 'GSM', 'SM', 'RTN', 'GSEQ', ...
  'HEE', 'HAE', 'HEEQ'};
for frameNo = 1 : length(frames)
  currFrame = frames{frameNo}
  % Convert to this frame
  cartTransPos = xcframe(time, cartPos     , [], [], 'GEO'    , currFrame);
  sphTransPos  = xcframe([]  , cartTransPos, [], [], currFrame, currFrame, ...
    'cart', 'sph' );
  
  % Display position in this frame
  disp(' ')
  disp(['Position in the ', currFrame, ' frame at ', datestr(time)])
  fprintf('Spherical coordinates: [%.2f, %.2f, %.2f] (deg, deg, km)\n', ...
    sphTransPos(1), sphTransPos(2), sphTransPos(3)*1e-3)
  fprintf('Cartesian coordinates: [%.2f, %.2f, %.2f] (km)\n', ...
    cartTransPos(1)*1e-3, cartTransPos(2)*1e-3, cartTransPos(3)*1e-3)
end
