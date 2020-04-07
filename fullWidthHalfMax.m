function [width,maxLocation,COM]=fullWidthHalfMax(x,y)

nInterp = 1000;
if length(x) < nInterp,
    xx = resample(x,nInterp,length(x));
    yy = resample(y,nInterp,length(y));
    x = xx;
    y = yy;
end;
maxY = max(y);
maxLocation = x(find(y==maxY));

halfstart = max(min(x),x(find(y(1:end-1) <= 0.5*maxY & y(2:end) > 0.5*maxY)));
halfend = min(max(x),x(find(y(1:end-1) >= 0.5*maxY & y(2:end) < 0.5*maxY)));

if isempty(halfstart), halfstart = NaN; halfend = NaN; end;
if isempty(halfend), halfstart = NaN; halfend = NaN; end;

if ~isnan(halfstart)
    halfstart = min([halfstart(end),max(halfstart(find(halfstart < maxLocation)))]);
    halfend = max([halfend(1),min(halfend(find(halfend > maxLocation)))]);
end
width = halfend - halfstart;
COM = sum(x.*y) ./ sum(y);
if isempty(width)
keyboard;
end;
end