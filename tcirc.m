function [p,c_rad]=tcirc(xydata,p_rad)
% [p,c_rad]=tcirc(xydata,p_rad) determines a p-value that a set
%   of Fourier components has a nonzero mean
%   and a confidence radius on the mean, via tcirc (Victor & Mast 1991)
%
% xydata: the data, xydata(1,:): real part, xydata(2,:), imaginary part
% p_rad: p-value for confidence radius; defaults to 0.05
%
% p: p-value for being significantly different from 0 (near 0: significant)
% c_rad: confidence radius (gets larger as p_rad gets smaller)
%
%  note that if p=p_rad, then this means borderline significance
%  and therefore c_rad^2=sum(xymean.^2)
%
%  J. Victor jdvicto@med.cornell.edu
%
% requires FINV, FCDF, in statistics toolbox
%
if (nargin<=1) p_rad=0.05; end
npts=size(xydata,2);
xymean=mean(xydata,2);
xycentered=xydata-repmat(xymean,1,npts);
xycov=(xycentered*xycentered')/(npts-1);
scov=xycov(1,1)+xycov(2,2);
%
fcirc=finv(1-p_rad,2,2*(npts-1));
c_rad=sqrt((1/npts)*fcirc*scov);
if (scov>0)
    p=1-fcdf(npts*sum(xymean.^2)/scov,2,2*(npts-1));
else
    if sum(xymean.^2)==0
        p=1; %zero mean, zero variance
    else
        p=0; %nonzero mean, zero variance
    end
end
return
