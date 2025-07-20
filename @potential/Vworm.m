function [V,F] = Vworm(cosAng,Vpm,varargin)
%--------------------------------------------------------------------------
        % Vworm performs the computation of a worm like bending force and
        % potential of polymerchain
        % input: 
        % cosAng - cosine of angles between subsequent subunits
        % Vpm - parameters, [kappa,kBT,Theta0], Theta0: designated angle
        % output:
        % V - potential value
        % F - currently not defined
        % optional:
        % see variable arguments
        %   See also VLennardJones,Vlinear
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % reference: https://doi.org/10.1021/acs.macromol.9b02437; modified
        % to any designated angle, not just pi
        % date: 2022/08/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('cosAng', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(cosAng,Vpm,varargin{:});
%--------------------------------------------------------------------------
%% computation
V=Vpm(1)*Vpm(2)*sqrt(abs(cos(Vpm(3))-cosAng));
F=nan;
%--------------------------------------------------------------------------
%figure;
% scatter3(r1(:,1),r1(:,2),r1(:,3),'filled');hold on;
% scatter3(r2(:,1),r2(:,2),r2(:,3),'filled');hold on;
% rescale=1000;
% quiver3(r1(:,1),r1(:,2),r1(:,3),F(:,1)*rescale,F(:,2)*rescale,F(:,3)*rescale,'linewidth',2);
end
