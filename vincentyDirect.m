%  [latitude2,longitude2,bearing2] = VINCENTYDIRECT(latitude1,longitude1,...
%     bearing1,distance,varargin)
%
%  DESCRIPTION
%  Estimates the latitude, longitude and final bearing of end point P2
%  (LATITUDE2,LONGITUDE2) from the geographic coordinates of the start point 
%  P1 (LATITUDE1,LONGITUDE1), its DISTANCE to the end point and the initial 
%  bearing BEARING1.
% 
%  This function uses Direct Vincenty formula, which bases its calculations 
%  on an elliptic model of the earth. Bearings are relative to north N and 
%  calculated considering travel path from P1 to P2 (P1->P2).
%
%  INPUT VARIABLES
%  - latitude1: Latitude of Point 1 [deg] 
%  - longitude1: Longitude of Point 1 [deg]
%  - distance: Distance between Point 1 and Point 2 [m]
%  - bearing1: initial bearing (direction P1->P2)
%  - ellipse (varargin{1}): character string to choose the Reference Ellipsoid. 
%    This parameter admits up to 25 different input strings between 
%    reference ellipsoids or datums.
%    
%    ellipse     DESCRIPTION
%    --------------------------------------------------------------------
%    'WGS84'     World Geodetic System 1984. Accurate general use (DEFAULT)
%    'WGS72'     World Geodetic System 1972. NASA, Dep. of Defense, marine
%                 seismic surveys.
%    'WGS66'     World Geodetic System 1960. General use
%    'WGS60'     World Geodetic System 1960. General use
%    'GRS80'     Geodetic Reference System 1980. North America and 
%                 Australia. Same precision as WGS84
%    'NAD83'     North American Datum 1983 (=GRS80). North America and 
%                 Australia
%    'GDA94'     Geocentric Datum of Australia 1994 (=GRS80). Australia 
%                 current Datum.
%    'CLK80'     Clarke 1880. Africa, France
%    'CLK66'     Clarke 1866. North America, Philippines
%    'NAD27'     North American Datum 1927 (=CLK66). North America
%    'AIR30'     Airy 1830. Great Britain
%    'MdAIR'     Modified Airy 1849. Ireland (Ireland 1965/75 Datum)
%    'AusNS'     Australian National Spheroid ANS. Australia
%    'AGD66'     Autralian Geodetic Datum 1966 (=AusNS). Autralia
%    'AGD84'     Autralian Geodetic Datum 1984 (=AusNS). Autralia 
%    'INTER'     International 1924 (a.k.a. Hayford 1909). Europe
%    'IAU68'     International Astronomical Union 1968
%    'IAU65'     International Astronomical Union 1965
%    'GRS67'     Geodetic Reference System 1967. South America
%    'MdGRS'     Modified Geodetic Reference System 1967. South
%                 America
%    'SAD69'     South American Datum 1969 (=MdGRS). South America
%    'KRASO'     Krasovsky 1940. East Europe / URSS
%    'ATS77'     Average Terrestrial System 1977. Canada Maritime Provinces
%    'EVRST'     Everest 1830. India, Burma, Pakistan, Afghanistan, Tailand
%    'BESSL'     Bessel 1841. Central Europe, Chile, Indonesia
%
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
%  - refEllip
%
%  FUNCTION CALLS
%  1) [latitude2,longitude2,bearing2] = vincentyDirect(latitude1,...
%       longitude1,bearing1,distance)
%     ¬ ellipse = 'WGS84', warn = 'on'
%  2) [latitude2,longitude2,bearing2] = vincentyDirect(latitude1,...
%       longitude1,bearing1,distance,ellipse)
%     ¬ warn = 'on'
%  3) [latitude2,longitude2,bearing2] = vincentyDirect(latitude1,...
%       longitude1,bearing1,distance,ellipse,warn)
%
%  CONSIDERATIONS & LIMITATIONS
%  - Vincenty formula is slower than Haversine but accurate at every 
%    distance (1 mm error withing considered elliptic earth model)
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong-vincenty.html
%  - Directorate of Overseas Surveys, "Survey Review", April 1975
%  - http://en.wikipedia.org/wiki/Vincenty's_formulae

%  VERSION 2.1
%  Date: 07 Jun 2015
%  Author: Guillermo Jimenez Arranz
%  - Added variable input argument 'warn' to disable the warning message
%    about coincident source and receiver points.
%
%  VERSION 2.0
%  Date: 02 dec 2014
%  Author: Guillermo Jimenez Arranz
%  - Enabled calculations for multiple positions (vector input). Improved 
%    performance against loops (for/while statements)
%  - Upgrade in comments
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  3 Jun 2014

function [latitude2,longitude2,bearing2] = vincentyDirect(latitude1,...
    longitude1,bearing1,distance,varargin)

% Variable Input Arguments
narginchk(4,6)
switch nargin
    case 4
        ellipse = 'WGS84';    
        warn = 'on';
    case 5
        ellipse = varargin{1};
        warn = 'on';   
    case 6
        ellipse = varargin{1};
        warn = varargin{2};
end

% Error Management (General)
if ~any(warn == [0 1])
    error('Variable input argument WARN must be [0 1] or logical')
end

% Manage Identical Points (P1=P2)
samePoint = (distance==0);
if any(samePoint) && warn
     warning(sprintf('One or more points share the same position\n'));  %#ok<SPWRN>
end

% General
latitude1 = latitude1*pi/180; % latitude Point 1 [rad]
longitude1 = longitude1*pi/180; % longitude Point 1 [rad]
bearing1 = bearing1*pi/180; % bearing Point 1 [rad]
[a,b,f] = refEllip(ellipse); % parameters for Reference Ellipsoid

% P2 Position and Final Bearing Calculation (latitude2,longitude2,bearing2)
tanU1 = (1-f)*tan(latitude1); 
cosU1 = cos(atan2(tanU1,1)); 
sigma1 = atan2(tanU1,cos(bearing1)); 
sinBeta = cosU1.*sin(bearing1);
cosSqBeta = (1-sinBeta).*(1+sinBeta); 
uSq = cosSqBeta*(a^2-b^2)/b^2; 
A = 1 + uSq.*(4096 + uSq.*(-768 + uSq.*(320 - 175*uSq)))/16384; 
B = uSq.*(256 + uSq.*(-128 + uSq.*(74 - 47*uSq)))/1024; 
sigma = distance./(b*A); % longitude in the auxiliary sphere

maxIte = 100; % maximum iteration
tol = 10^-12; % lambda tolerance (10^-12 implies 1 mm precision)
flag = 1; %loop access flag

while flag 
    o2SigmaM = 2*sigma1 + sigma; 
    sinSigma = sin(sigma); 
    cosSigma = cos(sigma); 
    cos2SigmaM = cos(o2SigmaM); 
    cosSq2SigmaM = cos2SigmaM.^2; 
    
    deltaSigma = B.*sinSigma.*(cos2SigmaM + 1/4*B.*(cosSigma.*(-1 ...
        + 2*cosSq2SigmaM) - 1/6*B.*cos2SigmaM.*(-3 + 4*sinSigma.^2).*(-3 ...
        + 4*cosSq2SigmaM)));
    sigmaP = sigma;
    sigma = distance./(b*A) + deltaSigma;
     
    err = abs(sigma-sigmaP); % approximation error
    maxIte = maxIte-1; % remaining iterations
    flag = (maxIte>0) && any(err>=tol); % loop access flag   
end

sinU1 = sin(atan2(tanU1,1));
cosSigma = cos(sigma);
sinSigma = sin(sigma);
cosBeta1 = cos(bearing1);
sinBeta1 = sin(bearing1);

lambda = atan2(sinSigma.*sinBeta1,cosU1.*cosSigma - sinU1.*sinSigma.*cosBeta1);
C = f/16*cosSqBeta.*(4+f*(4-3*cosSqBeta));
L = lambda - (1-C)*f.*sinBeta.*(sigma + C.*sinSigma.*(cos2SigmaM ...
    + C.*cosSigma.*(-1 + 2*cosSq2SigmaM)));

longitude2 = (L+longitude1)*180/pi; % longitude of point 2 [º]
latitude2 = atan2(sinU1.*cosSigma ...
    + cosU1.*sinSigma.*cosBeta1,(1-f)*sqrt(sinBeta.^2 ...
    + (sinU1.*sinSigma -  cosU1.*cosSigma.*cosBeta1).^2))*180/pi; % latitude of point 2 [º]
bearing2 = atan2(sinBeta,-sinU1.*sinSigma + cosU1.*cosSigma.*cosBeta1)*180/pi; % final bearing (P1-P2) [º]

