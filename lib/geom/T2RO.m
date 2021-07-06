function [R, O] = T2RO(T)
% Extracts R and O from T
% Andrea Ancillao
% 2020

R = T(1:3,1:3,:);
O = permute( T(1:3,4,:) , [3, 1, 2]); 
end
