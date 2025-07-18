function [obj] = Polymerization(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModActin'));
ip.addParameter('dt', [], @isnumeric); 
ip.addParameter('konG', [], @isnumeric);
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
%================================================================================
dt=ip.Results.dt;
if isempty(dt)
    dt=obj.pm.dt;
end
konG=ip.Results.konG;
if isempty(konG)
    konG=obj.pm.konG*ones(obj.var.Nbr,1);
end
for iBr=1:obj.var.Nbr
    if obj.var.capped(iBr)==false
    rdnNum=rand(1,1);
    if rdnNum<konG(iBr)*dt
        canPol=true;
    else
        canPol=false;
    end
    else
        canPol=false;
    end
    if canPol==true
        obj.var.n_coord=obj.var.n_coord+1;
        obj.var.lBr(iBr)=obj.var.lBr(iBr)+1;
        iSUinst=sum(obj.var.lBr(1:iBr,:)); %where new subunit is inserted
        iSUinBr_old=obj.var.iSUinBr;
        obj.var.iSUinBr=[iSUinBr_old(1:iSUinst-1);... %pre insert
                         iSUinBr_old(iSUinst-1);...
                         iSUinBr_old(iSUinst:end)]; %new subunit following its own branch
        coord_old=obj.var.coord;
        coord_new=ComMath.rotAxis([0;0;obj.pm.Dsu],obj.var.angBr(iBr,1),obj.var.angBr(iBr,2));
        coord_new=coord_old(iSUinst-1,:)+coord_new';
        obj.var.coord=[coord_old(1:iSUinst-1,:);...
                       coord_new;...
                       coord_old(iSUinst:end,:)];                       
    end
end
%================================================================================
end