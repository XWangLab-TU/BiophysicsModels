function [V] = ModClathrin_ModMembrane8V(mod,int_name,id_td_mem,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('int_name', @(x) ischar(x));
ip.addRequired('id_td_mem', @(x) islogical(x));
ip.addParameter('mex', [], @isobject);
ip.parse(mod,int_name,id_td_mem,varargin{:});
%----------------------------------------------------------------------------------------
%f.int_comp.(int_name)=zeros(3,3,n);
%--------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane, mod.i_mod.ModMemAdapter];  %1: clathrin, 2: membrane 3: MemAdapter
%--------------------------------------------------------------------------
%%
%%
[foot] = getFoot(mod.mod{i_mod(1)});
V=zeros(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(3)}.var.n_coord);
%V_rep=zeros(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(2)}.var.n_adp);
   for ic= 1:mod.mod{i_mod(1)}.var.n_coord
     %%
    for j= 1:mod.mod{i_mod(3)}.var.n_coord
        for k=1:3
            i=(ic-1)*3+k;
        if id_td_mem(i,j)==true
            V(i,j)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:),foot(k,:,ic),mod.TypForce.pm.(['k_' int_name]));
        end
        end
    end
   end
%fprintf('%.12f\n',sum(V));
end

function V=potential(r1,r2,k)
   V=0.5*k*(sum((r1-r2).^2,2));
end