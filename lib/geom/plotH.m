function varargout=plotH(scale,varargin)
% PLOTH Plot of a trasformation homogeneous matrix
%   h = PLOTH( scale, H01, H02, ... , H0n ) plot a single thetraedron associated to a
%   n homogeneous matrix H0i with i=1...n
%
%   h = PLOTH( scale, H01, H02, ... , H0n, ratio ) use a specified ratio
%       for the height of z axis, the default values is 3


if mod(length(varargin),2)
    ratio=varargin{end};
    varargin(end)=[];
    NRf=(nargin-2)/2;
else
    ratio=3;
    NRf=(nargin-1)/2;
end

H=cat(3,varargin{1:2:end});
c=varargin(2:2:end);
H=permute(H(1:3,:,:),[2,1,3]);
H(1:3,:,:)=H(1:3,:,:)*scale;
H(3,:,:)=H(3,:,:)/ratio;
H(1:3,:,:)=H(1:3,:,:)+H([4 4 4],:,:);
F=H([4,1,2,4,2,3,4,3,1],:,:);
F=reshape(F,[3,3,3,NRf]);

newplot;
view(3);
for i=1:NRf
    h(i)=patch(F(:,:,1,i),F(:,:,2,i),F(:,:,3,i),c{i});
end

switch nargout
    case 1
        varargout{1}=h;
    otherwise
end

    
