function [p,zhat1,zhat2,t2c] = T2circ2(z1,z2);
% function [p,zhat1,zhat2,t2c] = T2circ2(z1,z2);
%
% T2circ calculates the two-sample statistic T^2_circ from Victor and Mast 
% (1991). T^2_circ is designed to assess the whether two samples of mulitple 
% estimates of Fourier components are significantly different.  The 
% fundamental assumption of T^2_circ is that the real and imaginary
% components of the Fourier components are independent and
% Gaussian-distributed with equal variances.  See the original reference
% for more detail:
%
%  Victor, J.D., and Mast, J. (1991) A new statistic for steady-state evoked 
%     potentials. Electroenceph. Clin. Neurophysiol. 78, 378-388. 
%
% INPUTS:
%     z1, a vector or matrix of complex numbers (Fourier components). Each
%        row is assumed to be a separate estimate of the component.  If
%        multiple columns are present, t2c is calculated for each column
%        separately.
%     z2, a vector or matrix of complex numbers (Fourier components). Each
%        row is assumed to be a separate estimate of the component.  If
%        multiple columns are present, t2c is calculated for each column
%        separately. Note the number of columns of z2 must equal the number
%        of columns of z1.
%     alpha, a scalar, to determine the significance level for the
%        confidence limits. Defaults to 0.05.
%
% OUTPUTS:
%     p, a vector of probabilities that the population of estimates z1 could
%        come from a process with the same expected value as the population z2.  
%        If p<alpha, then the null hypothesis that E<z1> == E<z2> can be 
%        rejected for that probability level.
%     zhat1, the plug-in estimate of E<z1>, namely, mean(z1,1).
%     zhat2, the plug-in estimate of E<z2>, namely, mean(z2,1).
%     t2c, the T^2_circ statistic for each column of z; theoretically
%        distributed as F_2,_2*(# of estimates - 1).
%
% Function originally written by Andrew Hudson, 3/19/2005.


if nargin<2 
    error('T2circ2 requires at least two inputs.');
end

[M1,N1] = size(z1);
if M1 == 1 
    warning('T2circ2: you passed a row vector in; will assume that these are estimates of a single variable');
    z1 = z1';
    [M1,N1] = size(z1);
end
[M2,N2] = size(z2);
if M1 == 1 
    warning('T2circ2: you passed a row vector in; will assume that these are estimates of a single variable');
    z2 = z2';
    [M2,N2] = size(z2);
end
if (N1~=N2), error('T2circ2: z1 and z2 must have the same number of columns.'); end;

if (M1<3 | M2<3), error('T2circ2: Estimation of T^2_circ requires at least 3 estimates for both z1 and z2'); end;


zhat1 = mean(z1,1);
residuals1 = z1 - repmat(zhat1,[M1 1]);
sumSqResid1 = sum(abs(residuals1).^2,1);

zhat2 = mean(z2,1);
residuals2 = z2 - repmat(zhat2,[M2 1]);
sumSqResid2 = sum(abs(residuals2).^2,1);

t2c = (M1+M2-2)*abs(zhat1-zhat2).^2./(sumSqResid1+sumSqResid2);

p = 1 - fcdf(t2c*M1*M2/(M1+M2),2,2*M1+2*M2 - 4);

return;