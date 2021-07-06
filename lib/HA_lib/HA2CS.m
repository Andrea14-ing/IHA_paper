function T = HA2CS(n,s)
% Builds a CS based on an Helical Axis defined by n and s
% The procedure is similar to the ones described in (De Schutter, 2010) and (Mannel, 2004).
% T is the homogeneous matrix of the functional CS expressed in the same CS as n and s.
%
% Andrea Ancillao
% Leuven, 4/2021

nF = size(n,1);
[~, z] = mvarray(n); % ensure n is a unit vector
tmp = cross(z,s); % auxiliary direction

x = cross(n,tmp);
[~,x] = mvarray(x);
[~,y] = mvarray( cross(z,x));
T = nan(4,4,nF); % build matrix T
for k=1:nF
    O = s(k,:)';
    e1 = x(k,:)';
    e2 = y(k,:)';
    e3 = z(k,:)';
    T(:,:,k) = [ [e1,e2,e3,O] ; [0,0,0,1] ];
end

end