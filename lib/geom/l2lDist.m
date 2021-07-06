function [dist, c1, c2] = l2lDist(l1, l2)
% Calculate the distance between two lines along the common normal.
%
% [dist, c1, c2] = l2lDist(l1, l2)
%
% l1 and l2 are the two lines defined as [n, s] (direction and position)
% size [nF x 6]
% 
% dist is the distance along the common normal
% c1 point on line 1
% c2 point on line 2
%
% Andrea Ancillao
% Leuven 10-2020
% andrea.ancillao@hotmail.com


n1 = l1(:,1:3);
s1 = l1(:,4:6);
n2 = l2(:,1:3);
s2 = l2(:,4:6);

n = cross(n1,n2,2);
[~, n] = mvarray(n);

c1 = s1 + ( ( dot( (s2-s1) , cross(n2, n), 2) )./ dot(n1, cross(n2, n) , 2) ).*n1;
c2 = s2 + ( ( dot( (s1-s2) , cross(n1, n), 2) )./ dot(n2, cross(n1, n) , 2) ).*n2;

dist = p2pDist(c1,c2);
end


