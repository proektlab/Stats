function theta_bar = circularMean(theta);
% function theta_bar = circularMean(theta);
%
% Computes the mean direction from a vector of angles. Theta is assumed to
% be in radians.
%
if ~all(isreal(theta)), theta = angle(theta); end;
theta_bar = atan2(mean(sin(theta)),mean(cos(theta)));