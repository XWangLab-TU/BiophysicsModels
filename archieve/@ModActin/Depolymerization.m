function [obj] = Depolymerization(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModActin'));
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
%================================================================================
for iBr=1:obj.var.Nbr
    canDepol=true;
    if canDepol==true
        obj.var.n_coord=obj.var.n_coord-1;
        obj.var.lBr(iBr)=obj.var.lBr(iBr)-1;
        iSUinst=sum(obj.var.lBr(1:iBr,:)); %where new subunit is removed
        iSUinBr_old=obj.var.iSUinBr;
        obj.var.iSUinBr=[iSUinBr_old(1:iSUinst);... %pre delete
                         iSUinBr_old(iSUinst+2:end)]; %pro delete
        coord_old=obj.var.coord;
        obj.var.coord=[coord_old(1:iSUinst,:);...
                       coord_old(iSUinst+2:end,:)];                       
    end
end
%================================================================================
end