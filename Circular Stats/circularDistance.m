function circDist = circularDistance(X,Y)
% function circDist = circularDistance(X,Y);
%
% Computes the circular analog of Y - X.  X and Y can be complex numbers,
% in which case the phase angle of Y relative to X is returned, or they 
% can be angles in radians.
%
% INPUT:  X,Y - must match in dimensions, or one can be a scalar.
%
% OUTPUT: circDist, a (# of variables) element vector containing the angle
%               between variables in radians.
%
% Written by Andrew Hudson (c), 7/11/2006.
if ~all(isreal(X)), X = angle(X); end;
if ~all(isreal(Y)), Y = angle(Y); end;

circDist = atan2(sin(Y-X),cos(Y-X));