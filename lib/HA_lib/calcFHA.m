function [FHA, MFHA] = calcFHA(T2_1)
% Calculates: 
% 1) the Finite Helical Axis, FHA from the pose T of body 2 with respect to body 1. Implemented as described in (Spoor et al. 1980, Shekarforoush et al. 2018).
% 2) the Mean Finite Helical Axis, MFHA Axis, as described in (Ehrig et al. 2019).
% 3) a local coordinate system attached to the FHA and MFHA as described in (De Schutter, 2010).
%
% [FHA , MFHA] = calcFHA(T2_1)
%
% Inputs:
%       T2_1 = [4 x 4 x nF] Pose of body 2 with respect to body 1. Origin in [m]
%      
% Outputs:
%       FHA.n    [nFx3]  direction of the FHA in the same CS as the T2_1
%       FHA.s    [nFx3]  position of the FHA in the same CS as T2_1
%       FHA.T    [4 x 4 x nF] CS attached to the FHA
%       FHA.phi  [nFx1] angle of rotation about the FHA  (toghether with n, this is equivalent to the axis-angle representation of R)
%       FHA.t    [nFx1] displacement along the FHA
%
%       MFHA.n    [1x3]  direction of the MFHA in the same CS 
%       MFHA.s    [1x3]  position of the MFHA in the same CS
%       MFHA.T    [4 x 4] CS attached to the MFHA
%
% NOTE: The T is always built having the direction of the IHA as z-axis 
%
% Andrea Ancillao
% Leuven 04-2021
%
%
%
% References
%   C. W. Spoor and F. E. Veldpaus, “Rigid body motion calculated from spatial co-ordinates of markers,” J. Biomech., vol. 13, no. 4, pp. 391–393, Jan. 1980.
%	M. Shekarforoush, J. E. Beveridge, D. A. Hart, C. B. Frank, and N. G. Shrive, “Correlation between translational and rotational kinematic abnormalities and osteoarthritis-like damage in two in vivo sheep injury models,” J. Biomech., vol. 75, pp. 67–76, 2018.
% 	J. De Schutter, “Invariant Description of Rigid Body Motion Trajectories,” J. Mech. Robot., vol. 2, no. 1, p. 011004, 2010.
%   R. M. Ehrig and M. O. Heller, “On intrinsic equivalences of the finite helical axis, the instantaneous helical axis, and the SARA approach. A mathematical perspective,” J. Biomech., vol. 84, pp. 4–10, 2019
%
%
% DISCLAIMER:
% You running this script/function means you will not blame the author(s) 
% if this software breaks your stuff, destroys your computer or kills your cat. 
% This script/function is provided AS IT IS, without warranty of any kind. 
% Author(s) disclaim all implied warranties including, without limitation, 
% any implied warranties of merchantability or of fitness for a particular purpose. 
% The entire risk arising out of the use or performance of the sample scripts 
% and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever 
% (including, without limitation, damages for loss of business profits, business interruption, 
% loss of business information, or other pecuniary loss) arising out of the use of or inability to use 
% the script or documentation. Neither this script/function, nor any part of it may be republished 
% without author(s) express written permission. 
% Author(s) retain the right to alter this disclaimer at any time. 


%T1_2 = Hinv(T2_1); % to have the results expressed in {1}
T1_2 = T2_1;
nF = size(T1_2,3);


%% Calc FHA
% Procesure according to Spoor et al. 1980, Shekarforoush et al. 2018
threshold = 1/180*pi;
ntot = nan(nF,3);
stot = ntot;
ttot = nan(nF,1);
phitot = ttot;
validsamples = ones(nF,1);
for k=1:nF-1
    T1 = T1_2(:,:,k);
    T2 = T1_2(:,:,k+1);
    T = T1*Hinv(T2);
    [R, v] = T2RO(T);
    SinPhi = (1/2) * sqrt( (R(3,2)-R(2,3))^2 + (R(1,3)-R(3,1))^2 + (R(2,1)-R(1,2))^2 );
    CosPhi = (1/2) *( R(1,1) + R(2,2) + R(3,3) - 1);
    if SinPhi <= (sqrt(2)/2)
        phi = asin(SinPhi);
    else
        phi = acos(CosPhi);
    end
    if phi < threshold
        validsamples(k) = 0;
    end
    n = 1/(2*SinPhi) .* [R(3,2)-R(2,3) , R(1,3)-R(3,1) , R(2,1)-R(1,2)];
    t = n*v';
    s = cross( (-(1/2)*n) , cross(n , v) ) +  cross( (SinPhi/(2*(1-CosPhi))) * n , v);
    ntot(k,:) = n;
    stot(k,:) = s;
    ttot(k) = t;
    phitot(k) = phi;
end
% ntot(end+1,:) = ntot(end,:);
% stot(end+1,:) = stot(end,:);   
% ttot(end+1) = ttot(end);
% phitot(end+1) = phitot(end);


% LCS attached to the FHA 
T = HA2CS(ntot,stot);

% remove samples below threshold
% for k=1:nF
%   if ~validsamples(k)
%       ntot(k,:) = nan(1,3);
%       stot(k,:) = nan(1,3);
%       T(:,:,k) = nan(4);
%   end
% end

% threshold
thrsh = 0.5; %[deg]
for k=1:nF
    if (phitot(k).*180/pi) < thrsh
        ntot(k) = nan;
        stot(k) = nan;
        T(:,:,k) = nan(4,4);
    end
end

% store FHA data
FHA.n = ntot;
FHA.s = stot;
FHA.T = T;
FHA.phi = phitot;
FHA.t = ttot;



%% Calculate the MFHA (Ehrig 2019)
% The calculation is based on the all possible combinations (i,j) of FHAs in the
% set. 

% Search the optimal direction 
Rtot = zeros(3);
for i=1:nF-1 
    for j=(i+1):nF
        T1 = T1_2(:,:,i);
        T2 = T1_2(:,:,j);
        T = T1*Hinv(T2);
        [R3, v] = T2RO(T);
        Rtot = nansum( cat(3, Rtot , R3) , 3);
    end
end
[U,~,~] = svd(Rtot); 
Nopt = U(:,1)';


% Search the optimal pivot point
count = 0;
for i=1:nF-1
    for j=(i+1):nF
        count = count+1;
        T1 = T1_2(:,:,i);
        T2 = T1_2(:,:,j);
        T = T1 * Hinv(T2);
        [R, v] = T2RO(T);
        %         axang = rotm2axang(R3);
        %         theta = axang(:,3);
        %         weight = sin(theta/2)^2;
        %         k = axang(:,1:3); % kij direction
        
        SinPhi = (1/2) * sqrt( (R(3,2)-R(2,3))^2 + (R(1,3)-R(3,1))^2 + (R(2,1)-R(1,2))^2 );
        CosPhi = (1/2) *( R(1,1) + R(2,2) + R(3,3) - 1);
        if SinPhi <= (sqrt(2)/2)
            phi = asin(SinPhi);
        else
            phi = acos(CosPhi);
        end
        k = 1/(2*SinPhi) .* [R(3,2)-R(2,3) , R(1,3)-R(3,1) , R(2,1)-R(1,2)];
        c = cross( (-(1/2)*k) , cross(k , v) ) +  cross( (SinPhi/(2*(1-CosPhi))) * k , v);
        weight = sin(phi/2)^2;
        Qi(:,:,count) = weight.*(eye(3)-(k'*k));
        ci(count,:) = ( (weight.*c') )';
    end
end
Q = nansum(Qi,3);
Q2 = nansum(ci,1)';
Sopt = pinv(Q)*Q2;
Sopt = Sopt';




% LCS attached to the MFHA
Topt = HA2CS(Nopt,Sopt);


% store MFHA data
MFHA.n = Nopt;
MFHA.s = Sopt;
MFHA.T = Topt;


end














