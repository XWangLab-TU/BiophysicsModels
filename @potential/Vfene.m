function [V,F] = Vfene(r1,r2,Vpm,varargin)
%--------------------------------------------------------------------------
        % Vlinear performs the computation of FENE potential for simulating
        % polymer chains
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters, [k,R0]
        % output:
        % V - potential value
        % F - force on r1
        % optional:
        % see variable arguments
        %   See also VLennardJones,Vlinear
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % reference: https://doi.org/10.1016/j.procs.2017.05.152
        % date: 2022/08/03
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(r1,r2,Vpm,varargin{:});
%--------------------------------------------------------------------------
%% computation
r=sqrt(sum((r1-r2).^2,2));
V=-0.5*Vpm(1)*Vpm(2)^2*log(1-(r/Vpm(2)).^2);
F=Vpm(1)*r.*(r2-r1)./(1-(r/Vpm(2)).^2);
%--------------------------------------------------------------------------
%figure;
% scatter3(r1(:,1),r1(:,2),r1(:,3),'filled');hold on;
% scatter3(r2(:,1),r2(:,2),r2(:,3),'filled');hold on;
% rescale=1000;
% quiver3(r1(:,1),r1(:,2),r1(:,3),F(:,1)*rescale,F(:,2)*rescale,F(:,3)*rescale,'linewidth',2);
end
