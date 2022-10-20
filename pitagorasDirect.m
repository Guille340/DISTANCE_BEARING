%  [latitude2,longitude2,bearing2] = pitagorasDirect(latitude1,longitude1,...
%     bearing1,distance,varargin)
%
%  DESCRIPTION
%  Estimates the latitude, longitude and final bearing of end point P2
%  (LATITUDE2,LONGITUDE2) from the geographic coordinates of the start point 
%  P1 (LATITUDE1,LONGITUDE1), its DISTANCE to the end point and the initial 
%  bearing BEARING1.
%
%  This function approximates the area covered by two points in the earth
%  to a flat plane. This assumption simplifies the calculations (fast 
%  processing),but its accuracy is limited to short distances between 
%  points (<5 km). Bearings are relative to north N and calculated
%  considering travel path from point 1 to point 2 (P1->P2).
%
%  INPUT VARIABLES
%  - latitude1: Latitude of Point 1 [deg] 
%  - longitude1: Longitude of Point 1 [deg]
%  - distance: Distance between Point 1 and Point 2 [m]
%  - bearing1: initial bearing (direction P1->P2)
%  - warn (varargin{1}): TRUE for displaying a warning when the intput 
%    arguments include coincident start and end points.
%
%  OUTPUT VARIABLES
%  - latitude2: vector of latitudes of ending points (P2) [deg] 
%  - longitude2: vector of longitudes of ending points (P2) [deg]
%  - bearing2: vector of final bearings between start and end points (P1->P2).
%    BEARING2 = BEARING1 if DISTANCE = 0.
%
%  INTERNALLY CALLED FUNCTIONS
%  - below360
%
%  FUNCTION CALLS
%  1) [latitude2,longitude2,bearing2] = pitagorasDirect(latitude1,longitude1,bearing1,distance)
%     ¬ warn = 'on'
%  2) [latitude2,longitude2,bearing2] = pitagorasDirect(latitude1,longitude1,bearing1,distance,warn)
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong.html
%
%  See also PITAGORAS, BELOW360

%  VERSION 2.1
%  Date: 21 Dec 2016
%  Author: Guillermo Jimenez Arranz
%  - Added variable input argument 'warn' to disable the warning message 
%    about coincident source and receiver points.
%
%  VERSION 2.0
%  Date: 2 Dec 2014
%  Author: Guillermo Jimenez Arranz
%  - Changed the method used for bearing calculation, from Haversine to
%    Cartesian.
%  - Simplified calculation of latitude and longitude
%  - Enabled calculations for multiple positions (vector input). Improved 
%    performance against loops (for/while statements)
%  - Upgrade in comments
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  3 Jun 2014

function [latitude2,longitude2,bearing2] = pitagorasDirect(latitude1,longitude1,bearing1,distance,varargin)

% Variable Input Arguments
narginchk(4,5)
switch nargin
    case 4
        warn = 'on';
    case 5
        warn = varargin{1};
end

% Error Management (General)
if ~any(warn == [0 1])
    error('Variable input argument WARN must be [0 1] or logical')
end

% Manage Identical Points (P1=P2)
samePoint = (distance==0); % vector of logical values. Point1 and Point2 are different (0) or identical (1)
if any(samePoint) && warn
     warning(sprintf('One or more points share the same position\n'));  %#ok<SPWRN>
end

% General
R = 6378137; % earth radius [m]
latitude1 = latitude1*pi/180; % latitude of point 1 [rad]
longitude1 = longitude1*pi/180; % longitude of point 1 [rad]
bearing1 = bearing1*pi/180; % initial bearing (P1->P2) [rad]

% Geodetic Coordinates Calculation Point 2 (latitude2,longitude2)
deltaLat = distance/R.*cos(bearing1);
deltaLon = deltaLat.*tan(bearing1)./cos(latitude1 + deltaLat/2);
latitude2 = (latitude1 + deltaLat)*180/pi; % latitude point 2 [º]
longitude2 = (longitude1 + deltaLon)*180/pi; % longitude point 2 [º]

% Final Bearing Calculation Point 2 (bearing2)
bearing2 = below360(bearing1*180/pi,'deg'); % final bearing [deg]

end

