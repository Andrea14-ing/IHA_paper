function d = p2pDist(P1,P2)
% Distance between 2 points
% d = p2pDist(P1,P2)
%
% inputs: P1=(x,y,z) e p2=(x,y,z) [ nF x 3 ];
% or P1=(x,y) e p2=(x,y) [ nF x 2 ];

% Programma di 
% Andrea Ancillao 
% Roma, 23-12-2012
% Revision 
% Andrea Ancillao
% Leuven, 4/2021


% consistency check
if size(P1,2) ~= size(P2,2)
    error('Dimensions of inputs are not consistent');
end
if size(P1,1) ~= size(P2,1)
    error('Number of input frames is not consistent');
end

nD = size(P1,2);
nF = size(P1,1);
d = NaN (nF,1);

if nD == 3
    for k=1:nF
        x1 = P1(k,1);
        y1 = P1(k,2);
        z1 = P1(k,3);
        x2 = P2(k,1);
        y2 = P2(k,2);
        z2 = P2(k,3);
        d(k) = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 );
    end
elseif nD == 2
    for k=1:nF
        x1 = P1(k,1);
        y1 = P1(k,2);
        x2 = P2(k,1);
        y2 = P2(k,2);
        d(k) = sqrt( (x2-x1)^2 + (y2-y1)^2  );
    
    end
else
    error('Data dimension is wrong');
end
    
end


