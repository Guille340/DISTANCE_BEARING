%  [distance,bearing1,bearing2] = VINCENTY(latitude1,longitude1,latitude2,...
%     longitude2,varargin)
% 
%  DESCRIPTION
%  Estimates the DISTANCE, initial bearing BEARING1 and final bearing BEARING2
%  between two points in the earth (P1,P2) given by coordinates (LATITUDE1,
%  LONGITUDE1) and (LATITUDE2,LONGITUDE2).
%
%  This function uses Inverse Vincenty formula, which bases its calculations 
%  on an elliptic model of the earth. Bearings are relative to north N and 
%  calculated considering travel path from P1 to P2 (P1->P2). 
%
%  INPUT VARIABLES
%  - latitude1: vector of latitudes of starting points (P1) [deg] 
%  - longitude1: vector of longitudes of starting points (P1) [deg]
%  - latitude2: vector of latitudes of ending points (P2) [deg] 
%  - longitude2: vector of longitudes of ending points (P2) [deg]
%  - ellipse (varargin{1}): character string specifying the Reference Ellipsoid. 
%    This parameter admits any of the strings listed below. These strings
%    represent ellipsoids or datums.
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
%  - distance: Distance between Point 1 and Point 2 [m]
%  - bearing1: initial bearing (direction P1->P2). NaN if P1=P2
%  - bearing2: final bearing (direction P1->P2). NaN if P1=P2
%
%  INTERNALLY CALLED FUNCTIONS
%  - refEllip
%
%  FUNCTION CALLS
%  1) [distance,bearing1,bearing2] = vincenty(latitude1,longitude1,...
%       latitude2,longitude2)
%      ¬ ellipse='WGS84', warn = 'on'
%  2) [distance,bearing1,bearing2] = vincenty(latitude1,longitude1,...
%       latitude2,longitude2,ellipse)
%     ¬ warn = 'on'
%  3) [distance,bearing1,bearing2] = vincenty(latitude1,longitude1,...
%       latitude2,longitude2,ellipse,warn)
%
%  CONSIDERATIONS & LIMITATIONS
%  - Vincenty formula is slower than Haversine but accurate at every 
%    distance (1 mm error withing considered elliptic earth model)
%  - In revision 2.0 the function was modified to work with vectors. Due 
%    to the recursive nature of Vincenty approach, it is difficult to make 
%    an independent treatment of every value. The "while loop" ends either
%    when the number of iterations exceeds "maxIte" or when the error 
%    "err" for every pair of positions is lower than "tol". It means that 
%    the output parameters (distance,bearing1,bearing2) will be more accurate 
%    for some pair of positions than others. However, the error will always be 
%    <=tol.
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong-vincenty.html
%  - Directorate of Overseas Surveys, "Survey Review", April 1975
%  - http://en.wikipedia.org/wiki/Vincenty's_formulae
%
%  See also VINCENTYDIRECT, BELOW360

%  VERSION 2.1
%  Date: 7 Jun 2015
%  Author: Guillermo Jimenez Arranz
%  - Added variable input argument 'warn' to disable the warning message 
%    about coincident source and receiver points.
%
%  VERSION 2.0
%  Date: 02 dec 2014
%  Author: Guillermo Jimenez Arranz
%  - Enabled calculations for multiple positions (vector input). Improved
%    performance against loops (for/while statements). 
%  - Upgrade in comments
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  3 Jun 2014

function [distance,bearing1,bearing2] = vincenty(latitude1,longitude1,...
    latitude2,longitude2,varargin)

% Variable Input Arguments
narginchk(4,6)
switch nargin
    case 4
        ellipse = 'WGS84';
        warn = true;
    case 5
        ellipse = varargin{1};
        warn = true;
    case 6
        ellipse = varargin{1};     
        warn = varargin{2};
end
     
% Error Management (General)
if ~any(warn == [0 1])
    error('Variable input argument WARN must be [0 1] or logical')
end

% Manage Identical Points (P1=P2)
samePoint = (latitude1==latitude2) & (longitude1==longitude2); 
if any(samePoint) && warn
     warning('One or more points share the same position')
end

% General
latitude1 = latitude1*pi/180; % latitude Point 1 [rad]
longitude1 = longitude1*pi/180; % longitude Point 1 [rad] 
latitude2 = latitude2*pi/180; % latitude Point 2 [rad] 
longitude2 = longitude2*pi/180; % longitude Point 2 [rad] 
[a,b,f] = refEllip(ellipse); % parameters for Reference Ellipsoid

% Distance and Bearing Calculation (distance,bearing1,bearing2)
L = longitude2-longitude1; % difference of longitudes [m] 
U1 = atan((1-f)*tan(latitude1)); 
U2 = atan((1-f)*tan(latitude2)); 
sinU1 = sin(U1); 
cosU1=cos(U1); 
sinU2 = sin(U2); 
cosU2=cos(U2); 

lambda=L; % longitude in the auxiliary sphere
maxIte=100; % maximum iteration
tol=10^-12; % lambda tolerance (10^-12 implies 1 mm precision)
flag=1; % loop access flag

while flag
    sinLambda = sin(lambda); cosLambda=cos(lambda);
    sinSigma = sqrt((cosU2.*sinLambda).^2 ...
        + (cosU1.*sinU2-sinU1.*cosU2.*cosLambda).^2); 
    cosSigma = sinU1.*sinU2 + cosU1.*cosU2.*cosLambda; 
    sigma = atan2(sinSigma,cosSigma); 
    sinBeta = cos(U1).*cos(U2).*sinLambda./sinSigma; 
    cosSqBeta = 1-sinBeta.^2; 
    cos2SigmaM = cosSigma - 2*sinU1.*sinU2./cosSqBeta; 
    cos2SigmaM(isnan(cos2SigmaM)) = 0; 
    C = f/16*cosSqBeta.*(4+f*(4-3*cosSqBeta)); 
    
    lambdaP = lambda; 
    lambda = L + (1-C)*f.*sinBeta.*(sigma + C.*sinSigma.*(cos2SigmaM ...
        + C.*cosSigma.*(-1+2*cos2SigmaM.^2))); 
    
    err = abs(lambda-lambdaP); 
    maxIte = maxIte-1; 
    flag = (maxIte>0) && any(err>=tol); 
end

if maxIte == 0
    error('Iterative formula has failed to converge!')
end

uSq = cosSqBeta*(a^2-b^2)/b^2; 
A = 1 + uSq.*(4096 + uSq.*(-768 + uSq.*(320 - 175*uSq)))/16384;
B = uSq.*(256 + uSq.*(-128 + uSq.*(74 - 47*uSq)))/1024;
deltaSigma = B.*sinSigma.*(cos2SigmaM + 1/4*B.*(cosSigma.*(-1 ...
    + 2*cos2SigmaM.^2) - 1/6*B.*cos2SigmaM.*(-3 + 4*sinSigma.^2).*(-3 ...
    + 4*cos2SigmaM.^2)));

distance = b*A.*(sigma - deltaSigma); distance=round(distance*1000)/1000; % distance between points 1 and 2 [m] 
bearing1 = atan2(cosU2.*sinLambda,cosU1.*sinU2-sinU1.*cosU2.*cosLambda)*180/pi; % azimuth point 1 [º] 
bearing2 = atan2(cosU1.*sinLambda,-sinU1.*cosU2+cosU1.*sinU2.*cosLambda)*180/pi; % azimuth point 2 [º]

distance(samePoint) = 0; % if P1=P2 then distance=0
bearing1(samePoint) = NaN; % if P1=P2 then bearing1=NaN
bearing2(samePoint) = NaN; % if P1=P2 then bearing2=NaN
