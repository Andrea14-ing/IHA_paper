function h = plot_ellipse(c,v_a,v_b,nP,varargin)
% PLOTH Plot of an ellipse 
%   h = PLOT_ELLIPSE( c,a,b,nP ) plot nP points of an ellipse
t=[0:nP-1]'/(nP-1)*2*pi;
p=[cos(t) sin(t)];
p=p*([v_a(:),v_b(:)])';
p(:,1)=c(1)+p(:,1);
p(:,2)=c(2)+p(:,2);
if nargin>4
    h = plot(p(:,1),p(:,2),varargin{:});
else
    h = plot(p(:,1),p(:,2));
end
