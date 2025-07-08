function [id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepulsive(f,mod,int_name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('int_name', @(x) ischar(x));
ip.parse(f,mod,int_name,varargin{:}); 
%----------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane];  %1: clathrin, 2: membrane
%----------------------------------------------------------------------------------------
%% add repulsive force
id_feet=[49 50 51];
id_reps_m_c_leg=[];
id_all_c=1:52;
id_all_m=(1:mod.mod{i_mod(2)}.var.n_on_coord)';
id_leg=[[(1:8)';(25:32)'],[(9:16)';(33:40)'],[(17:24)';(41:48)']];
id_leg_feet=[id_leg;id_feet];
id_all_adp=(1:mod.mod{i_mod(2)}.var.n_adp)';
for ic= 1:mod.mod{i_mod(1)}.var.n_coord
    %%
    coord_cl=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(:,1:3),52,mod.mod{i_mod(1)}.var.O(:,:,ic))+mod.mod{i_mod(1)}.var.coord(ic,:);
    d_cl_to_mem=zeros(52,1);
    id_nearest_m=zeros(52,1);
    for isu=1:52
        [d_cl_to_mem(isu),id_nearest_m(isu)]=min(sum((coord_cl(isu,:)-mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_on_coord,:)).^2,2));
        %d_cl_to_mem(isu)=sqrt(d_cl_to_mem(isu));
    end
    id_nearest_m=mod.mod{i_mod(2)}.var.id_on_coord(id_nearest_m);
    for ileg=1:3
        [d_cl_mem_min,id_cl_to_mem]=min(d_cl_to_mem(id_leg_feet(:,ileg)));
        
        if (d_cl_mem_min<mod.mod{i_mod(2)}.pm.Vdw.rl_max)

                if id_leg_feet(id_cl_to_mem,ileg) <= 8 
                    id_leg_feet_tem=25;
                elseif id_leg_feet(id_cl_to_mem,ileg) <= 16
                    id_leg_feet_tem=33;
                elseif id_leg_feet(id_cl_to_mem,ileg) <= 24
                    id_leg_feet_tem=41;
                else
                    id_leg_feet_tem=id_leg_feet(id_cl_to_mem,ileg);
                end
                
                i_mem=id_nearest_m(id_leg_feet_tem);
                abc=[i_mem,mod.mod{i_mod(2)}.var.j_T(i_mem,1:2)];
                r_abc=mod.mod{i_mod(2)}.var.coord(abc,:);
                dir_abc=cross((r_abc(2,:)-r_abc(1,:)),(r_abc(3,:)-r_abc(1,:)));
                r_all=coord_cl(id_leg_feet(:,ileg),:)-mod.mod{i_mod(2)}.var.coord(i_mem,:);
                
                dot_all=zeros(17,1);
                for i=1:17
                dot_all(i)=dot(r_all(i,:),dir_abc);
                end
                
                if (numel(dot_all(dot_all<0))~=0) && (numel(dot_all(dot_all>0))~=0)
                    id_reps_m_c_leg=[id_reps_m_c_leg;id_nearest_m(id_leg_feet_tem),ic,ileg];
                end
                
        end
    end
end
end


