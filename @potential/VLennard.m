function[V,L]=VLennard(r1,r2,Epw,r0,varargin)
%==========================================================================
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('Epw', @(x) isnumeric(x));
ip.addRequired('r0', @(x) isnumeric(x));
ip.parse(r1,r2,Epw,r0,varargin{:});
%==========================================================================
%%computation
r=abs(r1-r2);
V=4*Epw*((r0./r)^12-(r0./r)^6);
%==========================================================================
end