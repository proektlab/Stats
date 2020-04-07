function [P] = twoSampleCircTest(X1,X2)
% function [P] = twoSampleCircTest(X1,X2)

for i=1:size(X1,2),
    x1 = X1(:,i);
    x2 = X2(:,i);
    
    nz1 = find(x1);
    nz2 = find(x2);
    
    x1 = x1(nz1)./abs(x1(nz1));
    x2 = x2(nz2)./abs(x2(nz2));
    
    n1 = size(x1,1);
    n2 = size(x2,1);
    n = n1 + n2;
    x = [x1;x2];
    fromFirst = 1:size(x1,1);
    fromSecond = size(x1,1)+(1:size(x2,1));
    [scores,idx] = sort(angle(x));
    uniformScores = 2*pi*(1:n)/n;

    if n1>n2
        beta = uniformScores(find(ismember(idx,fromFirst)));
        %R_sq = sum(cos(beta).^2) + sum(sin(beta).^2); % compute the resultant
        R_sq = sum(mean([cos(beta') sin(beta')]).^2);
        S_star = (1 - 1/(2*n1))*2*n1*R_sq + n1*R_sq.^2/2; % compute the corrected score
        nz = n1;
    else
        beta = uniformScores(find(ismember(idx,fromSecond)));
        %R_sq = sum(cos(beta).^2) + sum(sin(beta).^2); % compute the resultant
        R_sq = sum(mean([cos(beta') sin(beta')]).^2);
        S_star = (1 - 1/(2*n2))*2*n2*R_sq + n2*R_sq.^2/2; % compute the corrected score
        nz = n2;
    end
    P(i) = 1-chi2cdf(S_star,2);
end