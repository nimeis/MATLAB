function accels = gravAccel(pos1,pos2,mass2,G)

%gravVel(pos1,pos2,mass1,mass2,gravity)
%   Detailed explanation goes here
g = G; %6.67430 * 10^-11;
r = pdists(pos2,pos1);
unitVector = (pos2-pos1)./r;
accels = g.*mass2.*((1./r.^2).*(unitVector));

end

