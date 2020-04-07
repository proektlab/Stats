function [jkSE,jkHat] = standardErrorFromJackknifePopulation(jk);

N = size(jk,1);
jkHat = mean(jk,1);
jkSE = sqrt(((N-1)/N).*diag((jk - repmat(jkHat,[N 1 1]))'*(jk - repmat(jkHat,[N 1 1]))))';