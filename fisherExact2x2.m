function [p]=fisherExact2x2(table);
%function p=fisherExact(table);
%
%Gives you the p-value from a Fisher Exact test on a 2x2 table, arranged as a 2x2 matrix:
%			class 1		class 2
%	class 3		   a		   b
%	class 4		   c		   d
%
% This function uses the gammaln function to avoid having to calculate large factorials.
%
% DH 4/4/2002, based on inspiration from JDV
marginals=[sum(table,1)'; sum(table,2)];
if any(diag(table)==0)  % Swap it around to avoid weirdness on the diagonal...
   table=[table(2); table(4); table(1); table(3)];
end

table=table(:);
p=0; % gives you your p-value
N = sum(sum(table,1));
while (min(table(:))>=0)
   p=p+exp(sum([gammaln(marginals(:)+1); -gammaln(table(find(table))+1); -gammaln(N+1)]));
   table=table+[1 -1 -1 1]';
end
