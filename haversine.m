%  [distance,bearing1,bearing2] = HAVERSINE(latitude1,longitude1,latitude2,...
%      longitude2,varargin)
%
%  DESCRIPTION
%  Estimates the DISTANCE, initial bearing BEARING1 and final bearing BEARING2
%  between two points in the earth (P1,P2) given by coordinates (LATITUDE1,
%  LONGITUDE1) and (LATITUDE2,LONGITUDE2).
%
%  This function uses Haversine formula, which bases its calculations in a 
%  spherical model of the earth. Bearings are relative to north N and 
%  calculated considering travel path from P1 to P2 (P1->P2).
%
%  INPUT VARIABLES
%  - latitude1: vector of latitudes of starting points (P1) [deg] 
%  - longitude1: vector of longitudes of starting points (P1) [deg]
%  - latitude2: vector of latitudes of ending points (P2) [deg] 
%  - longitude2: vector of longitudes of ending points (P2) [deg]
%  - warn (varargin{1}): TRUE for displaying a warning when the intput 
%    arguments include coincident start and end points.
%
%  OUTPUT VARIABLES
%  - distance: vector of distances between start and end points (P1,P2) [m]
%  - bearing1: vector of initial bearings from start to end points (P1->P2). 
%    NaN if P1 = P2
%  - bearing2: vector of final bearings from start to end points (P1->P2). 
%    NaN if P1 = P2
%
%  INTERNALLY CALLED FUNCTIONS
%  - below360
%
%  FUNCTION CALLS
%  1. [distance,bearing1,bearing2] = haversine(latitude1,longitude1,...
%       latitude2,longitude2)
%     ¬ warn = 'on'
%  2. [distance,bearing1,bearing2] = haversine(latitude1,longitude1,...
%       latitude2,longitude2,warn)
%
%  CONSIDERATIONS & LIMITATIONS
%  - The Haversine formula is fast and accurate at long distances, but is 
%    less accurate than vincenty formula. The Closer the points to the equator
%    the greater the accuracy of this function.
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong.html
%  - http://mathforum.org/library/drmath/view/51879.html
%  - http://en.wikipedia.org/wiki/Haversine_formula
%
%  See also HAVERSINEDIRECT, BELOW360

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
%  - Upgrade in comments
%  - Modification of value returned by bearing1 and bearing2 when P1=P2. Changed
%    from an unknown value to NaN.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  3 Jun 2014

function [distance,bearing1,bearing2] = haversine(latitude1,longitude1,...
    latitude2,longitude2,varargin)

% Variable Input Arguments
narginchk(4,5)
switch nargin
    case 4
        warn = false;
    case 5
        warn = varargin{1};
end

% Error Management (General)
if ~any(warn == [0 1])
    error('Variable input argument WARN must be [0 1] or logical')
end

% Manage Identical Points (P1=P2)
samePoint = (latitude1==latitude2) & (longitude1==longitude2);
if any(samePoint) && warn
     warning('One or more points share the same position');
end

% General
latitude1 = latitude1*pi/180; % latitude Point 1 [rad]
longitude1 = longitude1*pi/180; % longitude Point 1 [rad]
latitude2 = latitude2*pi/180; % latitude Point 2 [rad]
longitude2 = longitude2*pi/180; % longitude Point 2 [rad]

% Distance Calculation (distance)
R = 6378137; % Earth radius [m]
deltaLat = latitude2-latitude1; % difference of latitudes [º]
deltaLon = longitude2-longitude1; % difference of longitudes [º]

h = sin(deltaLat/2).^2 + cos(latitude2).*cos(latitude1).*sin(deltaLon/2).^2; % function haversin(s/R)
distance = 2*R*asin(sqrt(h)); % distance between points 1 and 2

% Bearing Angle Calculation (bearing1,bearing2)
y = sin(deltaLon).*cos(latitude2);
x = cos(latitude1).*sin(latitude2) - sin(latitude1).*cos(latitude2).*cos(deltaLon);
bearing1 = atan2(y,x)*180/pi;
bearing1 = below360(bearing1,'deg');
bearing1(samePoint) = NaN;

y = sin(-deltaLon).*cos(latitude1);
x = cos(latitude2).*sin(latitude1) - sin(latitude2).*cos(latitude1).*cos(-deltaLon);
bearing2 = atan2(y,x)*180/pi + 180;
bearing2 = below360(bearing2,'deg');
bearing2(samePoint) = NaN;
