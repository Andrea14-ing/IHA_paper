function [m,v]=mvarray(a)
% Norm and unitvector associated to a vector [nF x 3]
% Output: m [nFx1], v [nFx3]
m = sqrt (nansum (a.^2,2));
nD = size(a,2);
v = a./repmat(m,[1,nD]);
end