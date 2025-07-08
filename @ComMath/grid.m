function [G] = grid(coord,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('coord', @(x) isnumeric(x));
ip.parse(coord,varargin{:});
%--------------------------------------------------------------------------
nDim=size(coord,2);
if (nDim>3) 
    error('supported up to 3D')
end
%%
if nDim==3
    x=coord(:,1);
    y=coord(:,2);
    z=coord(:,3);
elseif nDim==2
    x=coord(:,1);
    y=coord(:,2);
    z=0;
else
    x=coord(:,1);
    y=0;
    z=0;
end
sx=size(x);
sy=size(y);
sz=size(z);
G=zeros(sx(1)*sy(1)*sz(1),3);
for i=1:sx(1)
    G((i-1)*sy(1)*sz(1)+1:i*sy(1)*sz(1),1)=x(i);
end
for i=1:sx(1)
    for j=1:sy(1)
        G((i-1)*sy(1)*sz(1)+(j-1)*sz(1)+1:(i-1)*sy(1)*sz(1)+j*sz(1),2)=y(j);
    end
end
for i=1:sx(1)*sy(1)
    G((i-1)*sz(1)+1:i*sz(1),3)=z(1:sz(1));
end
end