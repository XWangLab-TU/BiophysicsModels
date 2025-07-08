function [V] = ModClathrin8V(mod,int_name, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('int_name', @(x) ischar(x));
ip.addParameter('mex', [], @isobject);
ip.addParameter('i_c_given', [], @isnumeric);
ip.parse(mod,int_name,varargin{:});
%----------------------------------------------------------------------------------------
%f.int_comp.(int_name)=zeros(3,3,n);
i_c_given=ip.Results.i_c_given;
%--------------------------------------------------------------------------
            i_mod=mod.i_mod.ModClathrin;  %1: clathrin
%%
[waist] = getWaist(mod.mod{i_mod});
[CtrlPt] = getCtrlPt(mod.mod{i_mod});
% u=zeros(mod.mod{i_mod}.var.n_coord,3);
% for i_c=1:mod.mod{i_mod}.var.n_coord
% u(i_c,:)=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c));
% end
V=zeros(mod.mod{i_mod}.var.n_coord,mod.mod{i_mod}.var.n_coord);
id_all=1:3;
%cos_psi_opt=cos(mod.mod{i_mod}.pm.psi);
if isempty(i_c_given)
for i_c1=1:mod.mod{i_mod}.var.n_coord-1
    for i_c2=i_c1+1:mod.mod{i_mod}.var.n_coord
        if ismember(i_c2,mod.mod{i_mod}.var.connect(:,1,i_c1))
            i_leg2=id_all(mod.mod{i_mod}.var.connect(:,1,i_c1)==i_c2);
            i_leg1=id_all(mod.mod{i_mod}.var.connect(:,1,i_c2)==i_c1);
            i_leg2=mod.mod{i_mod}.var.connect(i_leg2,2,i_c1);
            i_leg1=mod.mod{i_mod}.var.connect(i_leg1,2,i_c2);
            
            V_tem=potential(waist(i_leg1,:,i_c1),mod.mod{i_mod}.var.coord(i_c2,:),mod.TypForce.pm.(['k_' int_name])(1));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(waist(i_leg2,:,i_c2),mod.mod{i_mod}.var.coord(i_c1,:),mod.TypForce.pm.(['k_' int_name])(1));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(CtrlPt(i_leg1,:,i_c1),CtrlPt(i_leg2,:,i_c2),mod.TypForce.pm.(['k_' int_name])(2));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(CtrlPt(i_leg1+3,:,i_c1),CtrlPt(i_leg2+3,:,i_c2),mod.TypForce.pm.(['k_' int_name])(2));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
        end
    end
end
else
    %%
    i_c1=i_c_given;
    for i_c2=1:mod.mod{i_mod}.var.n_coord
        if (ismember(i_c2,mod.mod{i_mod}.var.connect(:,1,i_c1))) 
            i_leg2=id_all(mod.mod{i_mod}.var.connect(:,1,i_c1)==i_c2);
            i_leg1=id_all(mod.mod{i_mod}.var.connect(:,1,i_c2)==i_c1);
            i_leg2=mod.mod{i_mod}.var.connect(i_leg2,2,i_c1);
            i_leg1=mod.mod{i_mod}.var.connect(i_leg1,2,i_c2);
            V(i_c1,i_c2)=0;
            V(i_c2,i_c1)=0;
            V_tem=potential(waist(i_leg1,:,i_c1),mod.mod{i_mod}.var.coord(i_c2,:),mod.TypForce.pm.(['k_' int_name])(1));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(waist(i_leg2,:,i_c2),mod.mod{i_mod}.var.coord(i_c1,:),mod.TypForce.pm.(['k_' int_name])(1));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(CtrlPt(i_leg1,:,i_c1),CtrlPt(i_leg2,:,i_c2),mod.TypForce.pm.(['k_' int_name])(2));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
            %%
            V_tem=potential(CtrlPt(i_leg1+3,:,i_c1),CtrlPt(i_leg2+3,:,i_c2),mod.TypForce.pm.(['k_' int_name])(2));
            V(i_c1,i_c2)=V(i_c1,i_c2)+V_tem;
            V(i_c2,i_c1)=V(i_c2,i_c1)+V_tem;
        end
    end
end
%fprintf('%.12f\n',sum(V));
end

function V=potential(r1,r2,k)
   V=0.5*k*(sum((r1-r2).^2,2));
end