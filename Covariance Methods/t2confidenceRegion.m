function [sem,prediction,zscore,tsq] = t2confidenceRegion(data,p,isjk);

if nargin<3
   isjk=0;
end

[n,ax]=size(data);
if ax>2
   error(' Only 2-d data currently supported');
end

muHat = mean(data,1);
[pc,zscore,lambda,tsq]=princomp(data);
u=pc(:,1);
v=pc(:,2);
th = 0:pi/50:2*pi;
c = sqrt(ax*(n-1)/(n-ax))*sqrt(finv(p,ax,n-ax));
if ~isjk
   for t=1:length(th)
      prediction(t,:) = muHat + c*sqrt(lambda(1))*cos(th(t))*u'+ c*sqrt(lambda(2))*sin(th(t))*v';
      sem(t,:) = muHat + c*sqrt(lambda(1)/n)*cos(th(t))*u'+ c*sqrt(lambda(2)/n)*sin(th(t))*v';
   end
else
   for t=1:length(th)
      prediction(t,:) = muHat + c*sqrt(lambda(1))*cos(th(t))*u'+ c*sqrt(lambda(2))*sin(th(t))*v';
      sem(t,:) = muHat + c*sqrt((n-1)*lambda(1)/n)*cos(th(t))*u'+ c*sqrt((n-1)*lambda(2)/n)*sin(th(t))*v';
   end
end

