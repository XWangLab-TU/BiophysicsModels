function [obj] = SetVar(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isobject(x));
ip.addParameter('update', false, @islogical);
ip.addParameter('coordAssigned', [], @isnumeric);
ip.addParameter('case_init', 1, @isnumeric); 
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
update=ip.Results.update;
case_init=ip.Results.case_init;
%================================================================================
coordAssigned=ip.Results.coordAssigned;
% if isempty(coordAssigned)
    if case_init==1
        ctr=[0 0 0];
        r_rdn=1;
        ang_rdn=[1 2]*pi*1;
        coordSU=[];
        angBr=[];
        iSUinBr=[];
        for iBr=1:4
%             coordSU_ctr_tem=ctr+r_rdn*(rand(1,3)-0.5)*2;
%             coordSU_ctr_tem=ctr+[0 r_rdn*(rand(1,2)-0.5)*2];
            coordSU_ctr_tem=ctr+[0 cos(iBr*pi*0.5) sin(iBr*pi*0.5)];
            coordSU_tem=[0 0 obj.pm.Dsu];
%             ang_tem=ang_rdn.*(rand(1,2)-0.5)*2;
            ang_tem=[0.5*pi 1.5*pi];
            coordSU_tem = ComMath.rotAxis(coordSU_tem',ang_tem(1),ang_tem(2));
            coordSU=[coordSU; ...
                     coordSU_ctr_tem; ...
                     coordSU_tem'+coordSU_ctr_tem];
            angBr=[angBr; ang_tem];
            iSUinBr=[iSUinBr;iBr*ones(2,1)];
        end
    end
% else
%     coordSU=coordAssigned;
% end
Nsu=size(coordSU,1);
Nbr=size(angBr,1);
obj.var=struct('coord',coordSU,'n_coord',Nsu,'idMesh',[],...
               'iSUinBr',iSUinBr,'Nbr',Nbr,'lBr',zeros(Nbr,1),'angBr',angBr,'capped',false(Nbr,1),'iSUasArp',zeros(Nbr,1));
iSUinBr=1;
for iSU=1:Nsu
    if obj.var.iSUinBr(iSU)==iSUinBr
      obj.var.lBr(iSUinBr)=obj.var.lBr(iSUinBr)+1;
    else
      iSUinBr=iSUinBr+1;
      obj.var.lBr(iSUinBr)=obj.var.lBr(iSUinBr)+1;
    end
end
end