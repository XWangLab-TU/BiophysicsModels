function [G] = grid3(x,y,z,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x', @(x) isnumeric(x));
ip.addRequired('y', @(x) isnumeric(x));
ip.addRequired('z', @(x) isnumeric(x));
ip.parse(x,y,z,varargin{:});
%--------------------------------------------------------------------------
sx=size(x);
sy=size(y);
sz=size(z);
if (sx(2)~=1) || (sy(2)~=1) || (sz(2)~=1)
    error('input arrays must all be Nx1')
end
%%
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

