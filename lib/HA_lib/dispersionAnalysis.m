function out = dispersionAnalysis(ihas, refT, d, isplot)
% Calculate a reference plane and the intersection between the IHAs and
% this plane. 
% Calculate the confidence ellipse over the cloud of points on the plane 
% The referece plane is oriented having its normal as the direction of the
% AHA, and x,y axes corresponding to the x,y axes of the AHA CS.
%
% out = dispersionAnalysis(ihas, refT, d, isplot)
%
% Inputs:
%  
% ihas = [nF x 6] = pairs n,s defining the set of IHA axes
% refT  [4 x 4 x 1] is the CS associated to the reference axis. 
%   The direction of the reference axis is assumed to be the z-component of the R matrix 
% d [1] is the distance at which the plane should be calculated. Default = 1;
% isplot if =1 plots the ellipse. defauult = 0;
%
%
% Outputs:
% Data Structure:
%   The ellipse parameters
%       out.o center of the ellipse
%       out.va direction of the 1st axis
%       out.vb direction of the 2nd axis
%       out.ma magnitude of the 1st axis
%       out.mb magnitude of the 2nd axis
%
% The points on the plane expressed in CS plane:   
%       out.x
%       out.y
% 
% The CS of the plane as Tpl_0
%       out.Tpl 
%


% Andrea Ancillao, Ph.D.
% Leuven 10-2020
% Edit
% Andrea Ancillao, Ph.D.
% Leuven 07/2021


if nargin < 3
    d = 1;
    isplot = 0;
end
if nargin < 4
    isplot = 0;
   
end

% Get ASA CS info
[R,ASAo] = T2RO(refT);
ASAn = permute( R(1:3, 3,:) , [3, 1, 2]);  % take the z axis
ASAx = permute( R(1:3, 1,:) , [3, 1, 2]);
ASAy = permute( R(1:3, 2,:) , [3, 1, 2]);

% define plane
dist = d;    
Op = ASAo + dist.*ASAn;
Np = ASAn;

% CS of the plane
[~, e3] = mvarray(Np);  %z-axis oriented as the ASA n direction
e1 = ASAx; % x axis
e2 = ASAy;
Tp_1 = [[e1',e2',e3',Op'] ; [0,0,0,1]]; 


Nisa = ihas(:,1:3);
Sisa = ihas(:,4:6);
nF = size(ihas,1);

% the segment must be defined by 2 points thus:
P1 = Sisa - 100 * Nisa;
P2 = Sisa + 100 * Nisa;
      
Inters = nan(nF,3);
for k=1:nF
    if  ~isnan(ihas(k,1))
        [Inters(k,:),check] = plane_line_intersect(Np,Op,P1(k,:),P2(k,:));
        if check ~= 1
            disp(['Check axis plane intersection at frame: ', num2str(k), ' exit code: ',num2str(check)]);
        end
    end
end
               
    
T1_pl = Hinv(Tp_1);
PlaneInters_pl = transfPoint(Inters,T1_pl); % Tntercept seen in local plane

% calculate confidence ellipse
track = PlaneInters_pl(:,1:2); % remove firs dimension as all the points lie on the plane x-y
ref = nanmean(track);
for k = 1:size(track,1)
    if isnan(track(k,1))
        track(k,:) = ref;
    end
end
if ~isnan(track(1,1))
    [o, va, vb, ma, mb] = ellipse_fit(track,0.95);
else
    o = nan(1,2);
    va = o;
    vb = o;
    ma = nan;
    mb = nan;
end



out.o = o;
out.va = va;
out.vb = vb;
out.ma = ma;
out.mb = mb;

out.x = PlaneInters_pl(:,1);
out.y = PlaneInters_pl(:,2);

out.Tpl = Tp_1;


if isplot
   figure;
   hold on; grid on; box on; grid minor;
   nP = size(PlaneInters_pl,1);
   for k=1:nP
       plot(PlaneInters_pl(k,1), PlaneInters_pl(k,2) , 'b*');
   end
   plot(o(1),o(2) , 'r*', 'LineWidth',4);
   plot_ellipse(o,va*ma,vb*mb,1000,'r','LineWidth',2);
   title('Intersection with ref. plane');
   xlabel('x');
   ylabel('y');
   axis equal
end
     
end
