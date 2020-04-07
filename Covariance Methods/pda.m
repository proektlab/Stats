function [stats] = pda(traindata,trainclasses,testdata,testclasses,omega);
% function [stats] = pda(traindata,trainclasses,testdata,testclasses,omega);
%
% traindata, [Nsamples Mvariables]
% trainclasses, an Nsamples long vector
% omega, an [Mvariables Mvariables] penalty matrix

[Nsamples,Mvariables] = size(traindata);

if nargin<5 omega = []; end;
if prod(size(omega)) == 1
    lambda = omega;
    omega = [];
else
    lambda = [];
end;
if isempty(omega)
    if isempty(lambda), lambda = 0.1; end;
    temp = eye(Mvariables);
    temp = temp - circshift(temp,-1);
    Dk = temp(1:Mvariables-1,:);
    Dk_1 = Dk(1:size(Dk,1)-1,:);
    Dk_1 = Dk_1(:,1:size(Dk_1,2)-1);
    omega = lambda*Dk'*Dk_1'*Dk_1*Dk;
    omega = (omega+omega.')/2;
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


G = trainclasses(:);

allClasses = unique(trainclasses);
Jclasses = length(allClasses);

Y = zeros(Nsamples,Jclasses);
for i=1:Jclasses
    Y(find(G==allClasses(i)),i) = 1;
end;

grand_mean = mean(traindata,1);
H = traindata - repmat(grand_mean,[Nsamples 1]); % remove the mean response
sig_11 = (1/Nsamples)*Y'*Y; % diagonal matrix of class proportions (N_i/Ntot)
sig_22 = (1/Nsamples)*(H'*H+omega); % penalized covariance of the predictors
sig_22 = 0.5*(sig_22+sig_22.'); % minimize roundoff assymmetry
sig_12 = (1/Nsamples)*Y'*H; % covariation of predictors and class variables
p = diag(sig_11);


% we will solve the LDA problem by first doing a canonical correlation
% analysis...
[output,isig_11,isig_22] = ccaFromCovarianceMatrices(sig_11,sig_22,sig_12,Nsamples);

retain = max([max(find(output.P>.5)) 3]);
theta = output.a(:,1:retain);
B_cca = output.b(:,1:retain);
% B_cca = B_cca(:,1:Jclasses);
alphas = output.rho(1:retain);
B_os = B_cca*diag(alphas);
B_lda = B_cca*diag((1 - alphas.^2).^(-0.5));
M = isig_11*sig_12;

Hat = H*isig_22*H';
dof1 = trace(Hat);
dof2 = trace(Hat'*Hat);
dof = trace(Hat)-1;
residuals = (eye(Nsamples)-Hat)*Y*theta;
ASR = trace(residuals'*residuals)/Nsamples; % Average squared residuals
penalty = 2;
GCV = ASR/(1-(1+penalty*dof)/Nsamples)^2;
% AICc = Nsamples*log(sum(diag(residuals'*residuals))) + 2*dof1*Nsamples./(Nsamples - dof1 - 1);

train_projection = (H*B_lda).';
if ~isempty(testdata)
    test_projection = ((testdata-repmat(grand_mean,[size(testdata,1) 1]))*B_lda).';
else
    test_projection = [];
end
mean_projection = (M*B_lda).'; % column: class mean, row: lda dimension

p = diag(sig_11);
if ~isempty(testdata)
    for j=1:Jclasses
        train_distances(j,:) = sqrt(sum((train_projection - repmat(mean_projection(:,j),[1 size(train_projection,2)])).^2,1)) - 2*log(p(j));
        test_distances(j,:) = sqrt(sum((test_projection - repmat(mean_projection(:,j),[1 size(test_projection,2)])).^2,1)) - 2*log(p(j));
    end
else
    for j=1:Jclasses
        train_distances(j,:) = sqrt(sum((train_projection - repmat(mean_projection(:,j),[1 size(train_projection,2)])).^2,1)) - 2*log(p(j));
    end
    test_distances = [];
end;
[temp,train_class] = min(train_distances);
[temp,test_class] = min(test_distances);

train_true_id = length(find(trainclasses==allClasses(train_class)))/length(trainclasses);
train_false_id = 1 -  train_true_id;

test_true_id = length(find(testclasses==allClasses(test_class)))/length(testclasses);
test_false_id = 1 -  test_true_id;

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
stats.test_false_id = test_false_id;
stats.test_true_id = test_true_id;
stats.dof1 = dof1;
stats.dof2 = dof2;
stats.GCV = GCV;
stats.rho = alphas;
stats.r2classes = output.r2a;
stats.r2predictors = output.r2b;
stats.grand_mean = grand_mean;
stats.Pval = output.P;
