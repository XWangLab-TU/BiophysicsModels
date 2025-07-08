function [V,F] = Vlinear(r1,r2,Vpm,varargin)
%--------------------------------------------------------------------------
        % Vlinear performs the computation of a linear potential centered
        % at given location
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters, including amplitude and distance at potential minimum 
        % output:
        % V - potential value
        % F - force on r1
        % optional:
        % see variable arguments
        %   See also VinMembrane
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/03
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(r1,r2,Vpm,varargin{:});
%----------------------------------------------------------------------------------------
%% computation
r=sqrt(sum((r1-r2).^2,2));
if r==0
    u=rand(1,3);
else
    u=(r2-r1)/r;
end
V=Vpm(1)*abs(r-Vpm(2));
if r==Vpm(2)
    F=0;
elseif r>Vpm(2)
    F=Vpm(1)*u;
else
    F=-Vpm(1)*u;
end
%--------------------------------------------------------------------------    
end
