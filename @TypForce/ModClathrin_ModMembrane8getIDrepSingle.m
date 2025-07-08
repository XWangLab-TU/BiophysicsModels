function [id_reps_m_leg] = ModClathrin_ModMembrane8getIDrepSingle(f,mod,coord_new,O,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('coord_new', @(x) isnumeric(x));
ip.addRequired('O', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(f,mod,coord_new,O,varargin{:}); 
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane];  %1: clathrin, 2: membrane
%----------------------------------------------------------------------------------------
%% add repulsive force
id_feet=[49 50 51];
id_reps_m_leg=[];

id_leg=[[(1:8)';(25:32)'],[(9:16)';(33:40)'],[(17:24)';(41:48)']];
id_leg_feet=[id_leg;id_feet];

%%
    coord_cl=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(:,1:3),52,O)+coord_new;
%--------------------------------------------------------------------------
if ip.Results.plot_or_not==true
fig=figure('units','normalized','outerposition',[0 0 0.5 1]); 
x_lim=[-15 15];
plot(mod.mod{3},'linestyle','-','f',fig,'FaceAlpha',1);
% plot(mod.mod{1},'f',fig,'proxy_only',false,...
%      'simple',true,'iC_iLeg',[mod.mod{1}.var.n_coord 1;mod.mod{1}.var.n_coord 2;mod.mod{1}.var.n_coord 3]);
xlim(x_lim);ylim(x_lim);zlim(x_lim);
hold on;
scatter3(coord_cl(id_leg(:,1),1),coord_cl(id_leg(:,1),2),coord_cl(id_leg(:,1),3),30,'filled','markerfacecolor',[1 0 0]);hold on;
scatter3(coord_cl(id_leg(:,2),1),coord_cl(id_leg(:,2),2),coord_cl(id_leg(:,2),3),30,'filled','markerfacecolor',[0 1 0]);hold on;
scatter3(coord_cl(id_leg(:,3),1),coord_cl(id_leg(:,3),2),coord_cl(id_leg(:,3),3),30,'filled','markerfacecolor',[0 0 1]);hold on;
end
%--------------------------------------------------------------------------
    d_cl_to_mem=zeros(52,1);
    id_nearest_m=zeros(52,1);
    for isu=1:52
        [d_cl_to_mem(isu),id_nearest_m(isu)]=min(sqrt(sum((coord_cl(isu,:)-mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_on_coord,:)).^2,2)));
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
                    id_reps_m_leg=[id_reps_m_leg;id_nearest_m(id_leg_feet_tem),ileg];
                end
                
%                 if (id_leg(8,ileg)-id_cl_to_mem) < 0
%                     d_in_mem=id_leg(16,ileg)-id_cl_to_mem;
%                 else
%                     d_in_mem=id_leg(8,ileg)-id_cl_to_mem;
%                 end
%                 n_in_mem=[n_in_mem;d_in_mem];
        end
    end

end


