function h = plotPlane(n,O,scale,col,trasp)
% plots a plane as a square based on normal vector and origin
% n = [1x3]

%   $Revision: 1.0 $  $Date: 2019/05
%   $Author: A. Ancillao, Ph.D. 2019%
% 

if nargin < 5
    trasp = 0;
end


% v1 = rotx(90)*n'; % component on yz plane
% v2 = roty(90)*n'; % component on xz plane
% v3 = rotz(90)*n'; % component on xy plane
% 

% define a plane in 0
P1_0 = scale.*[1;1;0];
P2_0 = scale.*[1;-1;0];
P3_0 = scale.*[-1;-1;0];
P4_0 = scale.*[-1;1;0];
n_0 = [0;0;1];


CS1 = buildsdr([0,0,0],[0,0,1],[1,0,0]);
T1_0 = CS1;
T0_1 = Hinv(T1_0);

P1_1 = transfPoint(P1_0',T0_1)';
P2_1 = transfPoint(P2_0',T0_1)';
P3_1 = transfPoint(P3_0',T0_1)';
P4_1 = transfPoint(P4_0',T0_1)';


% P1 = O'+v1;
% P2 = O'-v1;

%% new sdr
[~ , i] = mvarray(n);
if i(1) ~= 1
    tmp = cross(i,[1,0,0]);
else
    tmp = cross(i,[0,1,0]);
end

[~ , k] = mvarray(tmp); 
[~ , j] = mvarray( cross( k , i ) ); 
R = cat(3,i,j,k); 
R = permute(R,[2,3,1]);
H = [ R , permute(O,[2,3,1]) ];
H(4,4,:) = 1; 
T2_0 = H;

P1_3 = transfPoint(P1_1',T2_0)';
P2_3 = transfPoint(P2_1',T2_0)';
P3_3 = transfPoint(P3_1',T2_0)';
P4_3 = transfPoint(P4_1',T2_0)';

corners = [P1_3' ; P2_3' ;P3_3' ;P4_3' ];
h = patch(corners(:,1) , corners(:,2), corners(:,3),col);
if trasp
    alpha(h,0.5);
end

end


function H = buildsdr(Centro,P1,P2)
centro2P1 = P1 - Centro;
centro2P2 = P2 - Centro;
[~ , i] = mvarray(centro2P1); 
tmp = cross(centro2P1,centro2P2);
[~ , k] = mvarray(tmp); 
[~ , j] = mvarray( cross( k , i ) );

R = cat(3,i,j,k); 
R = permute(R,[2,3,1]); 
H = [ R , permute(Centro,[2,3,1]) ];
H(4,4,:) = 1; 

end