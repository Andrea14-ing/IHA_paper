function v1 = transfVect(v0,R0_1)
%Rotation of a vector by matrix R
%   V0 vector in 0
%   R0_1 matrix form 0 to 1

%   $Date: 2019/02/01
%   $Author: A. Ancillao, 2019%


if size(v0,1) == 1 && size(R0_1,3) > 1
    nF = size(R0_1,3);
    v0 = repmat(v0,nF,1);
end
    
if size(R0_1,3) ==1 && size(v0,1) > 1
    nF = size(v0,1);
    R0_1 = repmat(R0_1,1,1,nF);
end

nF = size(v0,1);

if size(R0_1,3) == 1
    R0_1 = repmat(R0_1,1,1,nF);
end

v1 = nan(nF,3);
for k=1:nF
    v1(k,:) = ( R0_1(:,:,k) * v0(k,:)' )';   
end

