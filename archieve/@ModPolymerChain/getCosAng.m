function [obj] = getCosAng(obj,varargin)
%--------------------------------------------------------------------------
        % getCosAng performs the computation of Angle of @ModPolymerChain,
        % namely the angle between each neighboring pair of subunits
        % input: 
        % obj - a @ModPolymerChain object
        % optional:
        % see variable arguments
        %   See also SetParameter,SetVar
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/19
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModPolymerChain'));
ip.addParameter('idx', [], @isnumeric); %indicates only partial subunits for computation
ip.parse(obj,varargin{:});
%--------------------------------------------------------------------------
idx = ip.Results.idx;
update=false;
if isempty(idx)
    idx=1:obj.var.n_coord;
else
    update=true;
end
nIdx=numel(idx);
%--------------------------------------------------------------------------
Nmark=[0;cumsum(obj.var.nSubunit)]; %cumsum to get the last subunit number at the end of each chain
if update==false
    obj.var.cosAng=nan(obj.var.n_coord,1);
    obj.var.dir=nan(obj.var.n_coord,3);
    for iC=1:obj.var.nChain
        iS=Nmark(iC)+1:Nmark(iC+1)-2;
        a=obj.var.coord(iS+1,:)-obj.var.coord(iS,:);
        b=obj.var.coord(iS+2,:)-obj.var.coord(iS+1,:);
        aNorm=sqrt(sum(a.^2,2));
        bNorm=sqrt(sum(b.^2,2));
        obj.var.cosAng(iS+1)=sum(a.*b,2)./(aNorm.*bNorm);
        obj.var.dir(iS+1,:)=cross(-a./aNorm,b./bNorm,2);
        obj.var.dir(iS+1,:)=obj.var.dir(iS+1,:)./vecnorm(obj.var.dir(iS+1,:),2,2);
        %% test-------------------------------------------------------------
%         idR=1:18;
% %         idR=2;
%         R=obj.var.coord(2:end-1,:);
%         scatter3(R(idR,1),R(idR,2),R(idR,3));
%         hold on
%         F=-a;
%         xlim([-10 10]);ylim(xlim);zlim(xlim);
%         quiver3(R(idR,1),R(idR,2),R(idR,3),F(:,1),F(:,2),F(:,3),'linewidth',2);
%         F=b;
%         hold on
%         quiver3(R(idR,1),R(idR,2),R(idR,3),F(:,1),F(:,2),F(:,3),'linewidth',2);
%         F=cross(-a./aNorm,b./bNorm,2);
%         hold on
%         quiver3(R(idR,1),R(idR,2),R(idR,3),F(:,1),F(:,2),F(:,3),'linewidth',2);
%         [dot(F(1,:),a(1,:)),dot(F(2,:),a(2,:)),dot(F(3,:),a(3,:))]
%         [dot(F(1,:),b(1,:)),dot(F(2,:),b(2,:)),dot(F(3,:),b(3,:))]
    end
else
    idAvoid=[Nmark+1;Nmark];
    for iIdx=1:nIdx
        iS=idx(iIdx)-1;
        if ~ismember(iS+1,idAvoid)
            a=obj.var.coord(iS+1,:)-obj.var.coord(iS,:);
            b=obj.var.coord(iS+2,:)-obj.var.coord(iS+1,:);
            obj.var.cosAng(iS+1)=sum(a.*b,2)./(sqrt(sum(a.^2,2)).*sqrt(sum(b.^2,2)));
            obj.var.dir(iS+1,:)=cross(a,b);
            obj.var.dir(iS+1,:)=obj.var.dir(iS+1,:)./vecnorm(obj.var.dir(iS+1,:),2,2);
        end
    end
end
%--------------------------------------------------------------------------
end
%==========================================================================