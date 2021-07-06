function out = CorrectHAdir(HA,refdir)
% Compares the helical axis to a reference direction. If athe angle is
% greater than 90 deg, it inverts the HA.
% Use this to force the IHA to always point to the right of the subject.
% Use as reference any vector pointing roughly to the right of the subject.
%
% HA is a structure
%  .n [nF x 3]
%  .s [nF x 3]
%  .T [4 x 4 x nF]
%
% refdir is [nF x 6] or [nF x 3] or [1 x 6] or [1 x 3]
%
% Andrea Ancillao
% 29/04/2021
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


out = HA;
n = out.n;
N = refdir(:,1:3);

nF = size(n,1);
if size(N,1) == 1
    N = repmat(N,[nF,1]);
end


for k=1:nF
    angle  = vect2vectAng(n(k,:),N(k,:));
    if angle > 90
        out.n(k,:) = -out.n(k,:);
        % recalculate T
        ni = out.n(k,:);
        si = out.s(k,:);
        [~ , y] = mvarray( cross(ni,si));
        z = cross(ni,y);
        O = si';
        e1 = ni';
        e2 = y';
        e3 = z';
        out.T(:,:,k) = [ [e1,e2,e3,O] ; [0,0,0,1] ];
    end 
end

end


