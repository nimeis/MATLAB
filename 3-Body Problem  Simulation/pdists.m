function outputArg2 = pdists(pos,pos2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
deltas = pos2-pos;
deltas = deltas.^2;
outputArg2 = sum(deltas,2);
end