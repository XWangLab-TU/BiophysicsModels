function [S] = ellipsoid(a,c,varargin)
ip = inputParser;
ip.addRequired('a', @(x) isnumeric(x));
ip.addRequired('c', @(x) isnumeric(x));
ip.parse(a,c,varargin{:});
%==========================================================================
a = ip.Results.a;
c = ip.Results.c;
%==========================================================================
% a=b
%==========================================================================
idxAgtC=true(size(a));
idxAgtC(a<c)=false;
S=zeros(size(a));
e=zeros(size(a));
if sum(a==c)>0
    error('wrong, no sphere');
end

    e(~idxAgtC)=sqrt(1-a(~idxAgtC).^2./c(~idxAgtC).^2);
    S(~idxAgtC)=2*pi*a(~idxAgtC).^2.*(1+c(~idxAgtC)./a(~idxAgtC)./e(~idxAgtC).*asin(e(~idxAgtC)));

    e(idxAgtC)=sqrt(1-c(idxAgtC).^2./a(idxAgtC).^2);
    S(idxAgtC)=2*pi*a(idxAgtC).^2+pi*c(idxAgtC).^2./e(idxAgtC).*log((1+e(idxAgtC))./(1-e(idxAgtC)));

end