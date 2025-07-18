function [obj] = Severing(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModActin'));
ip.addParameter('dt', [], @isnumeric); 
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
%================================================================================
%%
dt=ip.Results.dt;
if isempty(dt)
    dt=obj.pm.dt;
end
NbrRem=0;
iBRrem=[];
iSUrem=[];
iSUall=(1:obj.var.n_coord)';
%%
for iBr=1:obj.var.Nbr
    rdnNum=rand(1,1);
    if rdnNum<obj.pm.ksev*dt
        canSev=true;
    else
        canSev=false;
    end
% for iBr=2:2
    if canSev==true
        NbrRem=NbrRem+1;
        idRem=obj.var.iSUinBr==iBr;
        nRem=numel(idRem(idRem==true));
        obj.var.n_coord=obj.var.n_coord-nRem;
        iBRrem=[iBRrem;iBr];
        iSUrem=[iSUrem;iSUall(idRem)];
    end
end
        obj.var.coord(iSUrem,:)=[];
        obj.var.iSUinBr(iSUrem)=[];
        obj.var.lBr(iBRrem)=[];
        obj.var.angBr(iBRrem,:)=[];
        obj.var.capped(iBRrem)=[];
%%        
obj.var.Nbr=obj.var.Nbr-NbrRem;
iBr=1;
if iBr<obj.var.Nbr
obj.var.iSUinBr(1:obj.var.lBr(iBr))=iBr;
idxStart=1;
idxEnd=obj.var.lBr(iBr);
for iBr=2:obj.var.Nbr
    idxStart=idxStart+obj.var.lBr(iBr-1);
    idxEnd=idxEnd+obj.var.lBr(iBr);
    obj.var.iSUinBr(idxStart:idxEnd)=iBr;
end
end
%================================================================================
end