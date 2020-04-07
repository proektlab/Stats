function omega = make2ndOrderSmoother(N,lambda);

if nargin<2 lambda = []; end;
if isempty(lambda), lambda = 0.1; end;

temp = eye(N);
temp = temp - circshift(temp,-1);
Dk = temp(1:N-1,:);
Dk_1 = Dk(1:size(Dk,1)-1,:);
Dk_1 = Dk_1(:,1:size(Dk_1,2)-1);
omega = lambda*Dk'*Dk_1'*Dk_1*Dk;