%  [latitude2,longitude2,bearing2] = VINCENTYDIRECT(latitude1,longitude1,...
%     bearing1,distance1,varargin)
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
%  - latitude1: latitude of start points P1 [deg] 
%  - longitude1: longitude of start points P1 [deg]
%  - bearing1: initial bearings (direction P1->P2)
%  - distance1: distances between start and end points (P1->P2) [m]
%
%  INPUT PROPERTIES
%  - warning: 0, 1 or logical. Use FALSE or 0 for omitting warnings.
%  - ellipsoid: character string representing the Reference Ellipsoid. 
%    This parameter admits up to 25 different input strings between 
%    reference ellipsoids or datums.
%    
%    ellipsoid   DESCRIPTION
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
%  OUTPUT VARIABLES
%  - latitude2: latitudes of end points (P2) [deg] 
%  - longitude2: longitudes of end points (P2) [deg]
%  - bearing2: final bearings between start and end points (P1->P2).
%    BEARING2 = BEARING1 if DISTANCE = 0.
%
%  INTERNALLY CALLED FUNCTIONS
%  - refEllip
%
%  FUNCTION CALLS
%  1. [latitude2,longitude2,bearing2] = vincentyDirect(latitude1,...
%       longitude1,bearing1,distance1)
%     * 'ellipsoid' = 'wgs84', 'warning' = TRUE
%  2. [latitude2,longitude2,bearing2] = vincentyDirect(...,<PROPERTY>,<VALUE>)
%
%  CONSIDERATIONS & LIMITATIONS
%  - Vincenty formula is slower than Haversine but accurate at every 
%    distance (1 mm error withing considered elliptic earth model)
%
%  REFERENCES
%  - http://www.movable-type.co.uk/scripts/latlong-vincenty.html
%  - Directorate of Overseas Surveys, "Survey Review", April 1975
%  - http://en.wikipedia.org/wiki/Vincenty's_formulae
%
%  See also vincenty, refEllip

%  VERSION 3.0
%  Date: 08 Apr 2023
%  Author: Guillermo Jimenez Arranz
%  Updates:
%  - Replaced variable input arguments WARN and ELLIPSOID with property/
%    value pairs 'warning' and 'ellipsoid'.
%  - Changed ELLIPSOID to ELLIPSE to avoid conflict with the function of
%    the same name.
%  - Changed DISTANCE to DISTANCE1 to avoid conflict with the function of
%    the same name.
%
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
    longitude1,bearing1,distance1,varargin)

    % Variable Input Arguments
    narginchk(4,8)
    nVarargin = nargin - 4;
    if rem(nVarargin,2)
        error('Property and value input arguments must come in pairs')
    end
    
    % Input Control
    [ellipse,warnFlag] = vincentyDirect_InputControl(latitude1,...
        longitude1,bearing1,distance1,varargin);

    % Manage Identical Points (P1=P2)
    samePoint = (distance1 == 0);
    if any(samePoint) && warnFlag
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
    sigma = distance1./(b*A); % longitude in the auxiliary sphere

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
        sigma = distance1./(b*A) + deltaSigma;

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
end

function [ellipse,warnFlag] = vincentyDirect_InputControl(latitude1,...
    longitude1,bearing1,distance1,varargin)
 
    % Initialise Default Parameters
    ellipse = 'wgs84';
    warnFlag = true;

    % Retrieve Input Variables
    varargin = varargin{1};
    nVarargin = length(varargin);
    for m = 1:2:nVarargin
        inputProperty = lower(varargin{m}); % case insensitive
        inputProperties = lower({'Ellipsoid','Warning'});
        if ~ismember(inputProperty,inputProperties)
            error('Invalid input property')
        else
            switch inputProperty
                case 'ellipsoid'
                    ellipse = varargin{m+1};
                case 'warning'
                    warnFlag = varargin{m+1};
            end
        end
    end
    
    % Error Control (general)
    if ~isnumeric(latitude1) || ~isvector(latitude1) ...
            || any(latitude1 < -90) || any(latitude1 > 90)
        error('LATITUDE1 must be a vector of numbers between -90 and 90.')
    end
    if ~isnumeric(longitude1) || ~isvector(longitude1) ...
            || any(longitude1 < -180) || any(longitude1 > 360)
        error('LONGITUDE1 must be a vector of numbers between -180 and 360.')
    end
    if ~isnumeric(bearing1) || ~isvector(bearing1) ...
            || any(bearing1 < -180) || any(bearing1 > 360)
        error('BEARING1 must be a vector of numbers between -180 and 360.')
    end
    if ~isnumeric(distance1) || ~isvector(distance1) || any(distance1 < 0)
        error('DISTANCE must be a vector of positive numbers.')
    end
    
    % Error Control (variable input arguments)
    ellipseIds = lower({'WGS84','WGS72','WGS66','WGS60','GRS80','NAD83',...
    'GDA94','AIR30','MdAIR','AusNS','AGD66','AGD84','INTER','IAU65','IAU68',...
    'GRS67','MdGRS','SAD69','CLK80','CLK66','NAD27','KRASO', 'ATS77',...
    'EVRST','BESSL'});
    if ~ischar(ellipse) || ~ismember(lower(ellipse),ellipseIds)
        ellipse = 'wgs84';
        warning(['ELLIPSOID is not a valid string identifier. '...
                'ELLIPSOID = ''WGS84'' will be used'])
    end
    if ~any(warnFlag == [0 1])
        warnFlag = true;
        warning(['WARNING must be [0 1] or logical. '...
            'Warnings will be displayed (WARNING = TRUE).'])
    end
end

