function [V,F] = VLennardJones(r1,r2,Vpm,varargin)
%--------------------------------------------------------------------------
        % Vlinear performs the computation of LennardJones between 2 points
        % in 3D
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters, [epsilon,sigma,order1], currently, the order of
        % the 2nd term is half of the order1
        % output:
        % V - potential value
        % F - force on r1
        % optional:
        % see variable arguments
        %   See also VinMembrane,Vfene,Vlinear
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
%--------------------------------------------------------------------------
%% computation
ord2=Vpm(3)/2;
r=sqrt(sum((r1-r2).^2,2));
V=4*Vpm(1)*((Vpm(2)./r).^Vpm(3)-(Vpm(2)./r).^ord2);
dVBYdr=4*Vpm(1)*(-Vpm(3)*Vpm(2)^Vpm(3)./r.^(Vpm(3)+1)+ord2*Vpm(2)^ord2./r.^(ord2+1)); 
F=(r2-r1)./r.*dVBYdr;
%--------------------------------------------------------------------------
%figure;
% scatter3(r1(:,1),r1(:,2),r1(:,3),'filled');hold on;
% scatter3(r2(:,1),r2(:,2),r2(:,3),'filled');hold on;
% rescale=1000;
% quiver3(r1(:,1),r1(:,2),r1(:,3),F(:,1)*rescale,F(:,2)*rescale,F(:,3)*rescale,'linewidth',2);
end
