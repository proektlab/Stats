function p = vonMisesPDF(x,mu,kappa)
% function p = vonMisesPDF(x,mu,kappa)
% Calculates the probability distribution function for a von Mises
% (circular normal) distributed variable with parameters mu and kappa.

p = exp(kappa*cos(x-mu))/(2*pi*besseli(0,kappa));

return;