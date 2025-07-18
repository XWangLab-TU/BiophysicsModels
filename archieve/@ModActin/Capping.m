function [obj] = Capping(obj, varargin)

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
for iBr=1:obj.var.Nbr
    rdnNum=rand(1,1);
    if rdnNum<obj.pm.kcap*dt
        canCap=true;
    else
        canCap=false;
    end
    if canCap==true
        obj.var.capped(iBr)=true;
    end
end
%================================================================================
end