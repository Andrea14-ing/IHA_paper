function T3 = Tcat(T1,T2)
% Concatenate homogeneus matrices T1 * T2 or orientation matrices R1 and R2 for all the frames
n1 = size(T1,3);
n2 = size(T2,3);
if n2 == 1 && n1 ~= 1
    T2 = repmat(T2,1,1,n1);
end

nF = size(T1,3);
nC = size(T1,2);
T3 = nan(nC,nC,nF);
for k=1:nF
    T3(:,:,k) = T1(:,:,k) * T2(:,:,k);
end
end