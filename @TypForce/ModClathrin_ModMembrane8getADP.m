function [id_td_mem] = ModClathrin_ModMembrane8getADP(mod,int_name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('int_name', @(x) ischar(x));
ip.addParameter('mex', [], @isobject);
ip.parse(mod,int_name,varargin{:});
%--------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane];  %1: clathrin, 2: membrane
%----------------------------------------------------------------------------------------
%f.int_comp.(int_name)=zeros(3,3,n);
%%
[foot] = getFoot(mod.mod{i_mod(1)});
if isempty(mod.TypForce.int_comp.(int_name))
    id_td_mem=false(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(2)}.var.n_adp);
else
    id_td_mem=false(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(2)}.var.n_adp);
    [nx,ny]=size(mod.TypForce.int_comp.(int_name){3});
    id_td_mem(1:nx,1:ny)=mod.TypForce.int_comp.(int_name){3};
end
%id_td_mem=false(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(2)}.var.n_adp);

% [id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepulsive(mod.TypForce,mod,i_mod,int_name,'ignore_adp',true);
% n_rep=size(id_reps_m_c_leg,1);

id_all=1:mod.mod{i_mod(2)}.var.n_adp;
for i= 1:mod.mod{i_mod(1)}.var.n_coord
%     if isempty(id_reps_m_c_leg)
%         is_memb=false;
%     else
%         [is_memb,id_mem]=ismember(i,id_reps_m_c_leg(:,2));
%     end
%     if is_memb==true  
%         [is_memb,id_adp]=ismember(id_reps_m_c_leg(id_mem,1),mod.mod{i_mod(2)}.var.id_adp);
%         if is_memb==true
%         id_td_mem((i-1)*3+id_reps_m_c_leg(id_mem,3),mod.mod{i_mod(2)}.var.id_adp(id_adp))=true;
%         end
%     else
    for k=1:3
        if mod.mod{i_mod(1)}.var.bound(i,k)==true
        %d_tem=sum(abs(id_td(k,:,i)-id_mem),2);
        d_tem=sum(abs(foot(k,:,i)-mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_adp,:)),2);
        [~,id_min]=min(d_tem);
%        if (d_min <=1)
%             [foot(k,:,i);
%              mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_adp(id_min),:)] 
%             pause;
            if sum(id_td_mem((i-1)*3+k,:))==0
            id_td_mem((i-1)*3+k,id_all(id_min))=true;
            end
%        end
        end
    end
%     end
end
end
