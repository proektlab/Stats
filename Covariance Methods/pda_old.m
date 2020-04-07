function [stats] = pda(traindata,classes,testdata,omega);
% function [stats] = pda(traindata,classes,testdata,omega);
%
% traindata, [Nsamples Mvariables]
% classes, an Nsamples long vector
% omega, an [Mvariables Mvariables] penalty matrix

[Nsamples,Mvariables] = size(traindata);

if nargin<4 omega = []; end;
if prod(size(omega)) == 1
    lambda = omega;
    omega = [];
else
    lambda = [];
end;
if isempty(omega)
    if isempty(lambda), lambda = 0.1; end;
    temp = eye(Mvariables); % diag(log(1+(1:Mvariables)));
    temp = temp - circshift(temp,-1);
    Dk = temp(1:Mvariables-1,:);
    Dk_1 = Dk(1:size(Dk,1)-1,:);
    Dk_1 = Dk_1(:,1:size(Dk_1,2)-1);
    omega = lambda*Dk'*Dk_1'*Dk_1*Dk;
    %omega = (omega+omega')/2;
end;

    
% check that omega meets its requirements
[on,om]=size(omega);
if (on~=Mvariables | om~=on)
    error('Omega must be a square matrix with as many columns as there are columns in the data matrix');
end

if any(any(omega ~= omega'))
    error('Omega must be symmetric');
end
e = svd(omega); 
if any(e<0&abs(e)>eps)
    keyboard;
    error('Omega must be nonnegative definite');
end;


G = classes(:);

allClasses = unique(classes);
Jclasses = length(allClasses);

Y = zeros(Nsamples,Jclasses);
for i=1:Jclasses
    Y(find(G==allClasses(i)),i) = 1;
end;

grand_mean = mean(traindata,1);
H = traindata - repmat(grand_mean,[Nsamples 1]); % remove the mean response
sig_11 = (1/Nsamples)*Y'*Y; % diagonal matrix of class proportions (N_i/Ntot)
sig_22 = (1/Nsamples)*(H'*H+omega); % penalized covariance of the predictors
sig_22 = 0.5*(sig_22+sig_22'); % minimize roundoff assymmetry
sig_12 = (1/Nsamples)*Y'*H; % covariation of predictors and class variables
sig_21 = sig_12';
p = diag(sig_11);

% we will solve the LDA problem by first doing a canonical correlation
% analysis...
[eig11,lambda11] = eig(sig_11); 
iss11 = eig11*inv(sqrt(lambda11))*inv(eig11);
[eig22,lambda22] = eig(sig_22); 
iss22 = eig22*inv(sqrt(lambda22))*inv(eig22);
iss22 = 0.5*(iss22+iss22'); % should be symmetric, but force down assymmetry due to numerical roundoff just in case

[t_1,Da,b_1] = svd(iss11*sig_12*iss22);

theta = iss11*t_1;
B_cca = iss22*b_1;
B_cca = B_cca(:,1:Jclasses);
alphas = diag(Da);
B_os = B_cca*diag(alphas);
B_lda = B_cca*diag((1 - alphas.^2).^(-0.5));

M = inv(sig_11)*sig_12;


S = H*eig22*inv(lambda22)*inv(eig22)*H'; % S is the so-called "hat" matrix
dof = trace(S)-1;
ASR = trace(eye(Jclasses-1)-Da(1:Jclasses-1,1:Jclasses-1).^2);
GCV = ASR/(1-(1+dof)/Nsamples)^2;

train_projection = B_lda'*H';
test_projection = B_lda'*(testdata-repmat(grand_mean,[size(testdata,1) 1]))';
mean_projection = B_lda'*M'; % column: class mean, row: lda dimension

p = diag(sig_11);
for j=1:Jclasses
    train_distances(j,:) = sqrt(sum((train_projection - repmat(mean_projection(:,j),[1 size(train_projection,2)])).^2,1)) - 2*log(p(j));
    test_distances(j,:) = sqrt(sum((test_projection - repmat(mean_projection(:,j),[1 size(test_projection,2)])).^2,1)) - 2*log(p(j));
end;
[temp,train_class] = min(train_distances);
[temp,test_class] = min(test_distances);

train_true_id = length(find(classes==allClasses(train_class)))/length(classes);
train_false_id = 1 -  train_true_id;

stats.B_lda = B_lda;
stats.B_os = B_os;
stats.B_cca = B_cca;
stats.theta = theta;
stats.train_projection = train_projection;
stats.test_projection = test_projection;
stats.mean_projection = mean_projection;
stats.train_id = allClasses(train_class);
stats.test_id = allClasses(test_class);
stats.train_false_id = train_false_id;
stats.train_true_id = train_true_id;
stats.dof = dof;
stats.GCV = GCV;
% stats.dof1 = dof1;
% stats.dof2 = dof2;
% stats.AIC = AIC;
% stats.AICc = AICc;
stats.grand_mean = grand_mean;
%stats.lambda = lambda;
for i=1:length(alphas)
    stats.Pval(i) = chi2cdf(-(Nsamples-1-0.5*(Jclasses+Mvariables+1)).*log(prod(1-alphas(1:i).^2)),(Jclasses-i+1)*(Mvariables-i+1));
end;
