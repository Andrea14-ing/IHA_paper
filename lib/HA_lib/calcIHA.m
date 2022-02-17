function [IHA , AHA] = calcIHA(tw)
% Calculates: 
% 1) the Instantaneous Helical Axis, IHA from the screw twist (tw) of body 2 with respect to body 1. Implemented as described in (Stokdijk et al 1999).
% 2) the Average Helical Axis, AHA base on the Optimal Pivot Point and Optimal direction Axis, as described in (Stokdijk et al 1999) but with weights on angular velocity as suggested by (Ehrig 2019).
% 3) a local coordinate system attached to the IHA and AHA. The calculation is similar to the ones described in (De Schutter, 2010) and (Mannel, 2004).
%
% [IHA , AHA] = calcIHA_Stokdijk(tw)
%
% Inputs:
%       tw = [nFx6] the screw twist of body 2 with respect to body 1. Units rad/s and m/s
%      
% Outputs:
%       IHA.n    [nFx3]  direction of the IHA in the same CS as the twist
%       IHA.s    [nFx3]  position of the IHA in the same CS as the twist [m]
%       IHA.T    [4 x 4 x nF] CS attached to the IHA
%
%       AHA.n    [1x3]  direction of the AHA in the same CS as the twist
%       AHA.s    [1x3]  position of the AHA in the same CS as the twist
%       AHA.T    [4 x 4] CS attached to the AHA
%
% NOTE: The T is always built having the direction of the IHA as z-axis 
%
%
% Andrea Ancillao
% Leuven 04-2021
%
%
% References
% 	M. Stokdijk, C. G. M. Meskers, H. E. J. Veeger, Y. A. De Boer, and P. M. Rozing, “Determination of the optimal elbow axis for evaluation of placement of prostheses,” Clin. Biomech., vol. 14, no. 3, pp. 177–184, 1999.
% 	J. De Schutter, “Invariant Description of Rigid Body Motion Trajectories,” J. Mech. Robot., vol. 2, no. 1, p. 011004, 2010.
% 	R. L. Lawrence, M. C. Ruder, R. Zauel, and M. J. Bey, “Instantaneous helical axis estimation of glenohumeral kinematics: The impact of rotator cuff pathology,” J. Biomech., vol. 109, p. 109924, 2020.
%   R. M. Ehrig and M. O. Heller, “On intrinsic equivalences of the finite helical axis, the instantaneous helical axis, and the SARA approach. A mathematical perspective,” J. Biomech., vol. 84, pp. 4–10, 2019
%   Mannel, H., Marin, F., Claes, L., & Dürselen, L. (2004). Establishment of a knee-joint coordinate system from helical axes analysis - A kinematic approach without anatomical referencing. IEEE Transactions on Biomedical Engineering, 51(8), 1341–1347. https://doi.org/10.1109/TBME.2004.828051
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



nF = size(tw,1);

% IHA
w = tw(:,1:3);
v = tw(:,4:6);
[wmod, n] = mvarray(w);
s = (cross(w,v)./(wmod.^2));


% define weights (based on w) (for optimal pivot point only)
wg = wmod./nanmax(wmod);
wg2 = LogisticFunc(wg,1,0.3,30); % logistic function

% apply the threshold and filter IHA
norig = n;
sorig = s;
trsh =  0.1 * nanmax(wmod); % rad/s
for k = 1:nF
    if wmod(k) <= trsh
        n(k,:) = nan(1,3);
        s(k,:) = nan(1,3);
    end
end

% LCS attached to the IHA
T = HA2CS(n,s);

% store IHA data
IHA.n = n;
IHA.s = s;
IHA.T = T;


% Optimal screw point - weighted
n = norig;
s = sorig;
Qi = nan(3,3,nF);
si = nan(nF,3);
for i=1:nF
       Qi(:,:,i) = wg2(i).*( eye(3)-(n(i,:)'*n(i,:)) );
       si(i,:) = ( Qi(:,:,i) * s(i,:)' )';
end
Q = nanmean(Qi,3);
Q2 = nanmean(si,1)';
Sopt = inv(Q)*Q2;
Sopt = Sopt';

% Optimal direction - weighted
ni = nan(3,3,nF);
for i=1:nF
       ni(:,:,i) = wg2(i).*(n(i,:)'*n(i,:) )';
end
ntot = nansum(ni,3);
[U,~,~] = svd(ntot); 
Nopt = U(:,1)';

% LCS attached to the AHA
Topt = HA2CS(Nopt,Sopt);


AHA.n = Nopt;
AHA.s = Sopt;
AHA.T = Topt;

end



