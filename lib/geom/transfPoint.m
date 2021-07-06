function v2 = transfPoint(v1,T1_2)
% Rotation and translation of a point by matrix T
%   V0 vector in 0
%   R0_1 matrix form 0 to 1

%   $Date: 2019/02/01
%   $Author: A. Ancillao, 2019%


nF = size(v1,1);

if size(T1_2,3) ==1
    T1_2 = repmat(T1_2,1,1,nF);
end

if size(v1,1) == 1 && size(T1_2,3) > 2
    nF = size(T1_2,3);
    v1 = repmat(v1,nF,1);
end
    

v2 = nan(nF,4);
for k=1:nF
    v2(k,:) = ( T1_2(:,:,k) * [v1(k,:) , 1]' )';   
end

v2(:,4) = [];

end
