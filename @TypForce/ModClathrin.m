function [f,V_tot,identifier] = ModClathrin(f,mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('idClathrinSub', [], @isnumeric);
ip.parse(f,mod,varargin{:});
%----------------------------------------------------------------------------------------
% if isempty(ip.Results.mex)
%     mex_avail=false;
% else
%     mex_avail=true;
%     mex=ip.Results.mex;
% end
mex_avail=true;
mex=mod.Mex;
int_name='ModClathrin';
%--------------------------------------------------------------------------
            i_mod=mod.i_mod.ModClathrin;  %1: clathrin
%--------------------------------------------------------------------------  
identifier=zeros(mod.n_mod,1);
identifier(mod.i_mod.ModClathrin)=1;
%--------------------------------------------------------------------------
%%
dr=0.000001;
if mex_avail==false
[V] = f.ModClathrin8V(mod,int_name);
V_tot=sum(sum(V));
%%
f.int_comp.(int_name)=cell(1,1);
f.int_comp.(int_name){1}=zeros(mod.mod{i_mod}.var.n_coord,6);
var_save=mod.mod{i_mod}.var;
for i_c= 1:mod.mod{i_mod}.var.n_coord
%--------------------------------------------------------------------------
        for k=1:3
%--------------------------------------------------------------------------   
           mod.mod{i_mod}.var.coord(i_c,k)=mod.mod{i_mod}.var.coord(i_c,k)+dr;
           [mod.mod{i_mod}.var] = setVar(mod.mod{i_mod},true);
           [V_alt] = f.ModClathrin8V(mod,int_name,'i_c_given',i_c);
           id_tem=V_alt~=0;
           f.int_comp.(int_name){1}(i_c,k)=-sum((V_alt(id_tem)-V(id_tem)))/dr;
           mod.mod{i_mod}.var=var_save;
%--------------------------------------------------------------------------
           mod.mod{i_mod}.var.a(i_c,k)=mod.mod{i_mod}.var.a(i_c,k)+dr;
           mod.mod{i_mod}.var.ang_a.Phi(i_c)=norm(mod.mod{i_mod}.var.a(i_c,:));
%            [mod.mod{i_mod}.var] = setVar(mod.mod{i_mod},true);
           var=mod.mod{i_mod}.var;
           var.O=mod.mod{i_mod}.Omega(var.a,var.ang_a.Phi);
           var.E=mod.mod{i_mod}.getE(var.a,var.ang_a.Phi);
           mod.mod{i_mod}.var=var;
           [V_alt] = f.ModClathrin8V(mod,int_name,'i_c_given',i_c);
           id_tem=V_alt~=0;
           f.int_comp.(int_name){1}(i_c,k+3)=-sum((V_alt(id_tem)-V(id_tem)))/dr;
           mod.mod{i_mod}.var=var_save;
%--------------------------------------------------------------------------
        end
end
else
    pmc=zeros(3,1);
    pmc(1) = mod.TypForce.pm.(['k_' int_name])(1);
    pmc(2) = mod.TypForce.pm.(['k_' int_name])(2);
    pmc(3) = dr;
    
    connect=[];
    for i=1:mod.mod{i_mod}.var.n_coord
        connect=[connect;mod.mod{1}.var.connect(:,:,i)];
    end
    
%%
    [f_c,V_tot]=mex.ModClathrin...
       (mod.mod{i_mod}.var.coord,...
        connect,...
        pmc,...
        mod.mod{i_mod}.var.coord_org,...
        mod.mod{i_mod}.var.a,...
        mod.mod{i_mod}.var.ctrl_pt_org);
%%
    f.int_comp.(int_name){1}=f_c;
end
f.int_V.ModClathrin=V_tot;
f.int_tot.ModClathrin=f.int_comp.ModClathrin;
end