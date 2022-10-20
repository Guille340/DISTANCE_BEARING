%  [a,b,f,ec,n] = REFELLIP(ellipsoid)
%
%  DESCRIPTION
%  Database of parameters for the most broadly used Reference Ellipsoids. 
%  A Reference Ellipsoid is a matemathical model of the earth that, along 
%  with its position relative to the center of the earth, establishes a 
%  Geographic Coordinate System or 'Datum'. Three parameters define the 
%  ellipsoid: radius at the equator A, radius at the poles B, and flattening
%  F. The same ellipsoid can be used for different Datums.
%
%  INPUT VARIABLES
%  - ellipsoid: character string specifying the Reference Ellipsoid. 25 types 
%    can be chosen using the following character strings.
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
%  OUTPUT VARIABLES
%  - a: radius at the equator for the selected ellipsoid [m]
%  - b: radius at the poles for the selected ellipsoid [m]
%  - f: flattening of the selected ellipsoid
%  - ec: eccentricity of the ellipsoid 
%  - n: third flattening of the ellipsoid
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The most accurate and widely used globally-applicable model for the 
%    earth ellipsoid is WGS-84, used in this script. Other ellipsoids 
%    offering a better fit to the local geoid include Airy (1830) in the UK, 
%    International 1924 in much of Europe, Clarke (1880) in Africa, GRS-67 
%    in South America, and many others. America (NAD83) and Australia (GDA) 
%    use GRS-80, functionally equivalent to the WGS-84 ellipsoid.
%
%  REFERENCES
%  - NIMA, "TR 8350.2, World Geodetic System 1984", January 2000, 3rd Ed.
%  - http://georepository.com/search/by-name/?query=&include_world=on
%
%  See also VINCENTY, VINCENTYDIRECT

%  VERSION 1.2
%  Date: 12 Jun 2014
%  Author: Guillermo Jimenez Arranz
%  - Added third flattening parameter as output (n)
%
%  VERSION 1.1
%  Date: 05 Jun 2014 
%  Author: Guillermo Jimenez Arranz
%  - Added new Ref. Ellipsoids and related Datums
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  03 Jun 2014

function [a,b,f,ec,n] = refEllip(ellipsoid)

switch ellipsoid
    case 'WGS84' % World Geodetic System 1984
        finv = 298.257223563; % inverse of ellipsoid flattening
        a = 6378137; % radius at equator [m]  
    case 'WGS72' % World Geodetic System 1972
        finv = 298.26; % inverse of ellipsoid flattening
        a = 6378135; % radius at equator [m]
    case 'WGS66' % World Geodetic System 1966
        finv = 298.25; % inverse of ellipsoid flattening
        a = 6378145; % radius at equator [m]
    case 'WGS60' % World Geodetic System 1960
        finv = 298.3; % inverse of ellipsoid flattening
        a = 6378165; % radius at equator [m]
    case {'GRS80','NAD83','GDA94'} % Geodetic Reference System 1980
        finv = 298.257222101; % inverse of ellipsoid flattening
        a = 6378137; % radius at equator [m]
    case 'AIR30' % Airy 1830
        finv = 299.3249646; % inverse of ellipsoid flattening
        a = 6377563.396; % radius at equator [m]
    case 'MdAIR' % Modified Airy 1849
        finv = 299.3249646; % inverse of ellipsoid flattening
        a = 6377340.189; % radius at equator [m]
    case {'AusNS','AGD66','AGD84'} % Autralian National Spheroid
        finv = 298.25; % inverse of ellipsoid flattening
        a = 6378160; % radius at equator [m]
    case 'INTER' % International 1924
        finv = 297; % inverse of ellipsoid flattening
        a = 6378388; % radius at equator [m]
    case 'IAU65' % International Astronomical Union 1965
        finv = 298.25; % inverse of ellipsoid flattening
        a = 6378160; % radius at equator [m]
    case 'IAU68' % International Astronomical Union 1968
        finv = 298.2472; % inverse of ellipsoid flattening
        a = 6378160; % radius at equator [m]
    case 'GRS67' % Geodetic Reference System 1967
        finv = 298.247167; % inverse of ellipsoid flattening
        a = 6378160; % radius at equator [m]
    case {'MdGRS','SAD69'} %  Modified GRS 1967
        finv = 298.25; % inverse of ellipsoid flattening
        a = 6378160; % radius at equator [m]
    case 'CLK80' % Clarke 1880
        finv = 293.465; % inverse of ellipsoid flattening
        a = 6378249.145; % radius at equator [m]
    case {'CLK66','NAD27'} % Clarke 1866
        finv = 294.9786982139; % inverse of ellipsoid flattening
        a = 6378206.4; % radius at equator [m]
    case 'KRASO' % Krasovsky 1940
        finv = 298.3; % inverse of ellipsoid flattening
        a = 6378245; % radius at equator [m]
    case 'ATS77' % Average Terrestrial System 1977
        finv = 298.257; % inverse of ellipsoid flattening
        a = 6378135; % radius at equator [m] 
    case 'EVRST' % Everest 1830
        finv = 300.8017; % inverse of ellipsoid flattening
        a = 6377276.345; % radius at equator [m] 
    case 'BESSL' % Bessel 1841
        finv = 299.1528128; % inverse of ellipsoid flattening
        a = 6377397.155; % radius at equator [m]        
    otherwise % World Geodetic System 1984 (DEFAULT MODEL)
        finv = 298.257223563; % inverse of ellipsoid flattening
        a = 6378137; % radius at equator [m]
end

f = 1/finv; % flattening of the ellipsoid
b = (1-f)*a; % radius at the poles [m]
ec = sqrt(1-(1-f)^2); % eccentricity of the ellipsoid
n = (a-b)/(a+b); % third flattening of the ellipsoid
