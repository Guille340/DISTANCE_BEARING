%  [latitude2,longitude2,bearing2] = HAVERSINEDIRECT(latitude1,longitude1,...
%     bearing1,distance,varargin)
% 
%  DESCRIPTION
%  Estimates the latitude, longitude and final bearing of end point P2
%  (LATITUDE2,LONGITUDE2) from the geographic coordinates of the start point 
%  P1 (LATITUDE1,LONGITUDE1), its DISTANCE to the end point and the initial 
%  bearing BEARING1.
%
%  This function uses the Haversine formula, which bases its calculations on 
%  a spherical model of the earth. Bearings are relative to north N and 
%  calculated considering a travel path from point 1 to point 2 (P1->P2).
%
%  INPUT VARIABLES
%  - latitude1: vector of latitudes of starting points (P1) [deg] 
%  - longitude1: vector of longitudes of starting points (P1) [deg] 
%  - distance: vector of distances between start and end points (P1->P2) [m]
%  - bearing1: vector of initial bearings between start and end points (P1->P2).
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
%  1. [latitude2,longitude2,bearing2] = haversineDirect(latitude1,...
%        longitude1,bearing1,distance)
%      ¬ warn = 'on'
%  2) [latitude2,longitude2,bearing2] = haversineDirect(latitude1,...
%        longitude1,bearing1,distance,warn)
%
%  CONSIDERATIONS & LIMITATIONS
%  - Haversine formula is fast and accurate at long distances, but is less 
%    accurate than vincenty formula. Closer the points to the equator
%    greater the accuracy of this function.
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong.html
%  - http://mathforum.org/library/drmath/view/51879.html
%  - http://en.wikipedia.org/wiki/Haversine_formula
%
%  See also HAVERSINE, BELOW360

%  VERSION 2.1
%  Date: 21 Dec 2016
%  Author: Guillermo Jimenez Arranz
%  - Added variable input argument 'warn' to disable the warning message 
%    about coincident source and receiver points.
%
%  VERSION 2.0
%  Date: 02 Dec 2014
%  Author: Guillermo Jimenez Arranz
%  - Enabled calculations for multiple positions (vector input). Improved 
%    performance against loops (for/while statements)
%  - Updated comments
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  3 Jun 2014

function [latitude2,longitude2,bearing2] = haversineDirect(latitude1,longitude1,bearing1,distance,varargin)

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

% Geodetic Coordinates Calculation P2 (latitude2,longitude2)
sinLat1 = sin(latitude1);
cosLat1 = cos(latitude1);
sinDR = sin(distance/R);
cosDR = cos(distance/R);

latitude2 = asin(sinLat1.*cosDR + cosLat1.*sinDR.*cos(bearing1));
y = sin(bearing1).*sinDR.*cosLat1;
x = cosDR - sinLat1.*sin(latitude2);
longitude2 = longitude1 + atan2(y,x);

% Final Bearing Calculation P2 (bearing2)
deltaLon = longitude1 - longitude2;
cosLat2 = cos(latitude2);
sinLat2 = sin(latitude2);
y = sin(deltaLon).*cosLat1;
x = cosLat2.*sinLat1 - sinLat2.*cosLat1.*cos(deltaLon);
if distance
    bearing2 = atan2(y,x)+pi; 
else
    bearing2 = bearing1;
end
latitude2 = latitude2*180/pi;
longitude2 = longitude2*180/pi;
bearing2 = below360(bearing2*180/pi,'deg');
