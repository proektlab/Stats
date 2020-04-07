function dataset = culloutliers(dataset,iqrmultiple)


if nargin<2, iqrmultiple = []; end;
if isempty(iqrmultiple), iqrmultiple = 1.5; end;
range = iqr(dataset);
m = median(dataset);
keep = find(dataset>m-iqrmultiple*range & dataset<m+iqrmultiple*range);
dataset = dataset(keep);