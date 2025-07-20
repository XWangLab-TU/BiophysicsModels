function [V,F] = ModPolymerChain(r1,r2,Vpm,varargin)
%--------------------------------------------------------------------------
        % ModPolymerChain performs the computation of potential and force 
        % between 2 subunits of ModPolymerChain in 3D
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters, including those for VLennardJones, Vfene and
        % Vcos
        % output:
        % V - potential value
        % F - force on r1
        % optional:
        % see variable arguments
        %   See also VinMembrane,Vfene,VLennardJones
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/18
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(r1,r2,Vpm,varargin{:});
%--------------------------------------------------------------------------
%% computation
[V,F]=potential.VLennardJones(r1,r2,Vpm(1:3));
[V1,F1]=potential.Vfene(r1,r2,Vpm(4:5));
V=V+V1;
F=F+F1;
% [V,F]=potential.Vfene(r1,r2,Vpm(4:5));
end
