function [ang, dist1, dist2] = compareAxes(ax1,ax2)
% Geometric comparison of two axes in Cartesian Space
% 
% [ang, dist1, dist2] = compareAxes(ax1,ax2)
%
% ax1 and ax2 axes defined as [n,s] = [nF x 6] direction and position
% note: ax1 OR ax2 may also be [1 x 6] (e.g. for the AHA). 
%
% Outputs:
% ang : angle between the two axes, in deg, range [0, 90]
% dist1: Linear distance between the two axes (along the common normal). Same unit as input s
% Dist2: Linear distance between the origins of the two axes
%
%
% Andrea Ancillao
% Leuven 10-2020
% Revision
% Andrea Ancillao
% Leuven 04/2021
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


nF1 = size(ax1,1);
nF2 = size(ax2,1);
if nF1 == 1 && nF2 > 1
    ax1 = repmat(ax1,[nF2,1]);
end
if nF2 == 1 && nF1 > 1
    ax2 = repmat(ax2,[nF1,1]);
end


%% angle between the two axes
n1 = ax1(:,1:3);
n2 = ax2(:,1:3);
ang = vect2vectAng( n1, n2);

% force angle to be [0, 90]
nF = size(ang,1);
for k=1:nF
    if ang(k) > 90
        ang(k) = 180 - ang(k);
    end
end


%% Distance along the common normal
dist1 = l2lDist(ax1, ax2);


%% Distance between origins
dist2 = p2pDist(ax1(:,4:6), ax2(:,4:6));

    
end
