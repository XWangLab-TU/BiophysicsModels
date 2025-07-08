function [V] = Vising(spin1,spin2,J,varargin)
%--------------------------------------------------------------------------
        
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('spin1', @(x) isnumeric(x));
ip.addRequired('spin2', @(x) isnumeric(x));
ip.addRequired('J', @(x) isnumeric(x));
ip.parse(spin1,spin2,J,varargin{:});
%----------------------------------------------------------------------------------------
%% computation
V=2*J*spin1.*spin2;
%--------------------------------------------------------------------------    
end
