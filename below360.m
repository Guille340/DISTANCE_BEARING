%  normAngle = BELOW360(angle,form)
%
%  DESCRIPTION
%  Calculates the corresponding angle between 0 and 360 degrees (or 0 and 2*pi 
%  radians) of an angle outside these limits (see also functions WRAPTO360, 
%  WRAPTO180, WRAPTO2PI and WRAPTOPI).
%
%  INPUT VARIABLES
%  - angle: angle to wrap. The units must be consistent with ?form
%  - form: string that describe the input units and type of wrapping
%    ¬ 'deg+': input in degrees, output in positive degrees (0 to 360)
%    ¬ 'deg': input in degrees, output in degrees (-180 to +180)
%    ¬ 'rad+': input in radians, output in positive radians (0 to 2*pi)
%    ¬ 'rad': input in radians, output in radians (-pi to +pi)
%
%  OUTPUT VARIABLES
%  - normAngle: wrapped angle [same units as ANGLE, same as FORM option]
%
%  INTERNALLY CALLED FUNCTIONS 
%  - None

%  VERSION 1.0: 09 Jun 2014
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com

function normAngle = below360(angle,form)

switch form
    case {'deg+' 'deg'}
        normAngle = rem(angle+360*ceil(abs(angle/360)),360);
        if strcmp(form,'deg')
            normAngle = normAngle - fix(normAngle/180)*360;
        end
    case {'rad+' 'rad'}
        normAngle = rem(angle+2*pi*ceil(abs(angle/(2*pi))),2*pi);
        if strcmp(form,'rad')
            normAngle = normAngle - fix(normAngle/pi)*2*pi;
        end      
end

