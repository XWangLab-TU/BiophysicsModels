function [V,F] = VwormDir(dir,Vpm,varargin)
%--------------------------------------------------------------------------
        % Vworm performs the computation of a worm like bending force and
        % potential of polymerchain
        % input: 
        % dir - direction of the subunits with respect to their neighbors
        % Vpm - parameters, only for magnitude
        % output:
        % V - potential value
        % F - currently not defined
        % optional:
        % see variable arguments
        %   See also VLennardJones,Vlinear, Vworm
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % to any designated angle, not just pi
        % date: 2022/08/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('dir', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(dir,Vpm,varargin{:});
%--------------------------------------------------------------------------
%% computation
% DirIdeal=0.3; %in radian
DirIdeal=cos(0.3); 
dotProd=sum(dir(1:end-1,:).*dir(2:end,:),2);
norm1=vecnorm(dir(1:end-1,:),2,2);
norm2=vecnorm(dir(2:end,:),2,2);
ang=acos(dotProd./norm1./norm2);
crossProd=cross(dir(1:end-1,:),dir(2:end,:),2);
hand=sum(crossProd.*crossProd(2,:),2);
idTem=hand<0;
ang(idTem)=-ang(idTem);
% ang(idTem)=2*pi-ang(idTem);
% n=numel(dirAlign)+1;
% Vdir1=sqrt(abs([dirAlign; nan]-DirIdeal)).*(1:n)';
% Vdir2=sqrt(abs([nan; dirAlign]-DirIdeal)).*(1:n)';
Vdir1=(abs([ang; nan]-DirIdeal)).^2;
Vdir2=(abs([nan; ang]-DirIdeal)).^2;
% Vdir1=sqrt(abs(cos([ang; nan])-DirIdeal));
% Vdir2=sqrt(abs(cos([nan; ang])-DirIdeal));
Vdir1(isnan(Vdir1))=0;
Vdir2(isnan(Vdir2))=0;
% Vdir=Vdir1+Vdir2;
Vdir=Vdir1;
V=Vpm(1)*sum(Vdir);
F=nan;
%--------------------------------------------------------------------------
%figure;
% scatter3(r1(:,1),r1(:,2),r1(:,3),'filled');hold on;
% scatter3(r2(:,1),r2(:,2),r2(:,3),'filled');hold on;
% rescale=1000;
% quiver3(r1(:,1),r1(:,2),r1(:,3),F(:,1)*rescale,F(:,2)*rescale,F(:,3)*rescale,'linewidth',2);
end
