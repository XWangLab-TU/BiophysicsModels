function [f,V_tot] = ModClathrin_ModMemAdapter_ModMembrane(f,mod,int_name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('int_name', @(x) ischar(x));
% ip.addRequired('A', @(x) isa(x,'ModMembrane'));
% ip.addRequired('B', @(x) isobject);
ip.addParameter('mex', [], @isobject);
ip.addParameter('update', true, @islogical);
ip.parse(f,mod,int_name,varargin{:}); 
%----------------------------------------------------------------------------------------
update=ip.Results.update;
if isempty(ip.Results.mex)
    mex_avail=false;
else
    mex_avail=true;
    mex=ip.Results.mex;
end
%----------------------------------------------------------------------------------------
%%
i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane, mod.i_mod.ModMemAdapter];  %1: clathrin, 2: membrane, 3: ModMemAdapter
%----------------------------------------------------------------------------------------
%%
if update==true
    f.int_comp.(int_name)=cell(4,1);
    f.int_comp.(int_name){1}=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
    f.int_comp.(int_name){2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);
    f.int_comp.(int_name){3}=false(mod.mod{i_mod(1)}.var.n_coord*3,mod.mod{i_mod(3)}.var.n_coord);
    id_all_adp=(1:mod.mod{i_mod(3)}.var.n_coord)';
    id_tem=mod.mod{i_mod(3)}.var.id_ModClathrin(:,1)>0;
    id_adp=id_all_adp(id_tem);
    id_c=mod.mod{i_mod(3)}.var.id_ModClathrin(id_tem,1);
    id_leg=mod.mod{i_mod(3)}.var.id_ModClathrin(id_tem,2);
    n_adp_c=numel(id_c);
    for i_adp=1:n_adp_c
        f.int_comp.(int_name){3}((id_c(i_adp)-1)*3+id_leg(i_adp),id_adp(i_adp))=true;
    end
    f.int_comp.(int_name){4}=zeros(mod.mod{3}.var.n_coord,1);
else
    f.int_comp.(int_name){1}=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
    f.int_comp.(int_name){2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);
    f.int_comp.(int_name){4}=zeros(mod.mod{3}.var.n_coord,1);
end
%%
dr=0.000001;
if mex_avail==false
[V] = mod.TypForce.ModClathrin_ModMemAdapter_ModMembrane8V(mod,int_name,f.int_comp.(int_name){3});
[foot] = getFoot(mod.mod{i_mod(1)});
mod_save=mod;
id_feet=[49 50 51];
%%
for i= 1:mod.mod{i_mod(1)}.var.n_coord
%--------------------------------------------------------------------------
    for j= 1:mod.mod{i_mod(3)}.var.n_coord
        for k=1:3
            id_tem=(i-1)*3+k;
        if (f.int_comp.(int_name){3}(id_tem,j)==true) %&& (mod.mod{i_mod(1)}.var.bound(i,k)==true)
%--------------------------------------------------------------------------   
           V_alt=V;
           for m=1:3
           mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),m)=mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),m)+dr;
           %[V_alt] = mod.TypForce.ModClathrin_ModMembrane8V(mod,i_mod,int_name,f.int_comp.(int_name){3});
           V_alt(id_tem,j)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:),foot(k,:,i),mod.TypForce.pm.(['k_' int_name]));  
           f_tem=-(V_alt(id_tem,j)-V(id_tem,j))/dr;
           %f_tem=f_tem/(mod.mod{i_mod(2)}.var.val(mod.mod{i_mod(3)}.var.id_ModMembrane(j))+1);
           id_jT=mod.mod{i_mod(2)}.var.j_T(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:);
           id_jT=id_jT(~isnan(id_jT));
           f.int_comp.(int_name){2}(mod.mod{i_mod(3)}.var.id_ModMembrane(j),m)=f.int_comp.(int_name){2}(mod.mod{i_mod(3)}.var.id_ModMembrane(j),m)+f_tem;
           f.int_comp.(int_name){2}(id_jT,m)=f.int_comp.(int_name){2}(id_jT,m)+f_tem;
           V_alt(id_tem,j)=V(id_tem,j);
           mod=mod_save;
           end
%--------------------------------------------------------------------------
           for m=1:3
%--------------------------------------------------------------------------   
           mod.mod{i_mod(1)}.var.coord(i,m)=mod.mod{i_mod(1)}.var.coord(i,m)+dr;
           [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);           
           foot_tem=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet(k),1:3),1,mod.mod{i_mod(1)}.var.O(:,:,i))+mod.mod{i_mod(1)}.var.coord(i,:);
           %[V_alt1] = mod.TypForce.ModClathrin_ModMembrane8V(mod,i_mod,int_name,f.int_comp.(int_name){3});
           V_alt(id_tem,j)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:),foot_tem,mod.TypForce.pm.(['k_' int_name]));
           f.int_comp.(int_name){1}(i,m)=f.int_comp.(int_name){1}(i,m)-((V_alt(id_tem,j))-(V(id_tem,j)))/dr;
           V_alt(id_tem,j)=V(id_tem,j);
           mod=mod_save;
%--------------------------------------------------------------------------
           mod.mod{i_mod(1)}.var.a(i,m)=mod.mod{i_mod(1)}.var.a(i,m)+dr;
           mod.mod{i_mod(1)}.var.ang_a.Phi(i)=norm(mod.mod{i_mod(1)}.var.a(i,:));
           [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);
           foot_tem=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet(k),1:3),1,mod.mod{i_mod(1)}.var.O(:,:,i))+mod.mod{i_mod(1)}.var.coord(i,:);
           %[V_alt1] = mod.TypForce.ModClathrin_ModMembrane8V(mod,i_mod,int_name,f.int_comp.(int_name){3});
           V_alt(id_tem,j)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:),foot_tem,mod.TypForce.pm.(['k_' int_name]));
           f.int_comp.(int_name){1}(i,m+3)=f.int_comp.(int_name){1}(i,m+3)-((V_alt(id_tem,j))-(V(id_tem,j)))/dr;
           V_alt(id_tem,j)=V(id_tem,j);
           mod=mod_save;
%--------------------------------------------------------------------------
           end
%--------------------------------------------------------------------------
%fprintf('%d,%d\n',i,j);pause;

           dV_adp=inf(mod.mod{i_mod(2)}.pm.n_val_max,1);
           i_m=mod.mod{i_mod(3)}.var.id_ModMembrane(j);
           j_T=mod.mod{i_mod(2)}.var.j_T(i_m,:);
           for m=1:mod.mod{i_mod(2)}.pm.n_val_max
               if ~isnan(j_T(m))
                   mod.mod{i_mod(3)}.var.id_ModMembrane(j)=j_T(m);

                   dV_adp(m)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),:),foot(k,:,i),mod.TypForce.pm.(['k_' int_name])); 

                   dV_adp(m)=dV_adp(m)-V(id_tem,j);
                   mod=mod_save;
               end
           end

%           hold on;scatter3(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),1),mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),2),mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(3)}.var.id_ModMembrane(j),3),30,'filled','markerfacecolor',[0 1 1]);
           %%
           [dV_adp_min,j_T_min] = min(dV_adp);
           if dV_adp_min<0
               f.int_comp.(int_name){4}(mod.mod{i_mod(3)}.var.id_ModMembrane(j))=j_T(j_T_min);
           end
%--------------------------------------------------------------------------    
        end
%--------------------------------------------------------------------------        
        end
%--------------------------------------------------------------------------      
    end
%--------------------------------------------------------------------------
end
V_tot=sum(sum(V));
%==========================================================================
%% add repulsive force
[id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepulsive(f,mod,int_name);
n_reps=size(id_reps_m_c_leg,1);
V=zeros(n_reps,1);
for i_reps=1:n_reps
    im=id_reps_m_c_leg(i_reps,1);
    ic=id_reps_m_c_leg(i_reps,2);
    ileg=id_reps_m_c_leg(i_reps,3);
    V(i_reps)=potential(mod.mod{i_mod(2)}.var.coord(im,:),foot(ileg,:,ic),mod.TypForce.pm.(['k_' int_name]));
end
V_alt=V;
%--------------------------------------------------------------------------
for i_reps=1:n_reps
    im=id_reps_m_c_leg(i_reps,1);
    ic=id_reps_m_c_leg(i_reps,2);
    ileg=id_reps_m_c_leg(i_reps,3);
%--------------------------------------------------------------------------
    for m=1:3
        mod.mod{i_mod(2)}.var.coord(im,m)=mod.mod{i_mod(2)}.var.coord(im,m)+dr;
        V_alt(i_reps)=potential(mod.mod{i_mod(2)}.var.coord(im,:),foot(ileg,:,ic),mod.TypForce.pm.(['k_' int_name]));
        f_tem=-(V_alt(i_reps)-V(i_reps))/dr;
        %f_tem=f_tem/(mod.mod{i_mod(2)}.var.val(im)+1);
        id_jT=mod.mod{i_mod(2)}.var.j_T(im,:);
        id_jT=id_jT(~isnan(id_jT));
        f.int_comp.(int_name){2}(im,m)=f.int_comp.(int_name){2}(im,m)+f_tem;
        f.int_comp.(int_name){2}(id_jT,m)=f.int_comp.(int_name){2}(id_jT,m)+f_tem;
        V_alt(i_reps)=V(i_reps);
        mod=mod_save;
    end
%--------------------------------------------------------------------------
%%
    for m=1:3     
        mod.mod{i_mod(1)}.var.coord(ic,m)=mod.mod{i_mod(1)}.var.coord(ic,m)+dr;
        [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);
        foot_tem=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet(ileg),1:3),1,mod.mod{i_mod(1)}.var.O(:,:,ic))+mod.mod{i_mod(1)}.var.coord(ic,:);
        V_alt(i_reps)=potential(mod.mod{i_mod(2)}.var.coord(im,:),foot_tem,mod.TypForce.pm.(['k_' int_name]));
        f.int_comp.(int_name){1}(ic,m)=f.int_comp.(int_name){1}(ic,m)-(V_alt(i_reps)-V(i_reps))/dr;
        V_alt(i_reps)=V(i_reps);
        mod=mod_save;
        %--------------------------------------------------------------------------
        mod.mod{i_mod(1)}.var.a(ic,m)=mod.mod{i_mod(1)}.var.a(ic,m)+dr;
        mod.mod{i_mod(1)}.var.ang_a.Phi(ic)=norm(mod.mod{i_mod(1)}.var.a(ic,:));
        [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);
        foot_tem=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet(ileg),1:3),1,mod.mod{i_mod(1)}.var.O(:,:,ic))+mod.mod{i_mod(1)}.var.coord(ic,:);
        V_alt(i_reps)=potential(mod.mod{i_mod(2)}.var.coord(im,:),foot_tem,mod.TypForce.pm.(['k_' int_name]));
        f.int_comp.(int_name){1}(ic,m+3)=f.int_comp.(int_name){1}(ic,m+3)-(V_alt(i_reps)-V(i_reps))/dr;
        V_alt(i_reps)=V(i_reps);
        mod=mod_save;
        %--------------------------------------------------------------------------
    end
%--------------------------------------------------------------------------
end
V_tot=V_tot+sum(V);
%--------------------------------------------------------------------------
else
    pmc=zeros(2,1);
    pmc(1) = mod.TypForce.pm.(['k_' int_name]);
    pmc(2) = dr;
    
    [id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepulsive(f,mod,int_name);
    
    j_T = mod.mod{i_mod(2)}.var.j_T; j_T(isnan(j_T)) = 0;
%%
    [f_c,f_m,V_tot,f_adp]=mex.ModClathrin_ModMemAdapter_ModMembrane...
       (mod.mod{i_mod(1)}.var.coord,...
        mod.mod{i_mod(2)}.var.coord,...
        mod.mod{i_mod(3)}.var.id_ModMembrane,...
        id_reps_m_c_leg,...
        j_T,...
        mod.mod{i_mod(2)}.var.id_on_coord,...
        pmc,...
        floor(f.int_comp.(int_name){3}+0.5),...
        mod.mod{i_mod(2)}.var.n_node',...
        mod.mod{i_mod(1)}.var.coord_org,...
        mod.mod{i_mod(1)}.var.a);
%%
    f.int_comp.(int_name){1}=f_c;
    f.int_comp.(int_name){2}=f_m;
    f.int_comp.(int_name){4}=floor(f_adp+0.5)+1;
end
%%
plot_or_not=false;
if plot_or_not==true
fig=figure;
plot(mod.mod{2},'linestyle',':','f',fig,'FaceAlpha',1);
col=[zeros(mod.mod{1}.var.n_coord,1) ones(mod.mod{1}.var.n_coord,1) ones(mod.mod{1}.var.n_coord,1)];
%col=rand(mod.mod{1}.var.n_coord,3);
plot(mod.mod{1},'f',fig,'proxy_only',false,'col',col,...
     'simple',true,'iC_iLeg',[8 2;7 2]);
for i_reps=1:n_reps
    foot_reps=foot(id_reps_m_c_leg(i_reps,3),:,id_reps_m_c_leg(i_reps,2));
    m_reps=mod.mod{i_mod(2)}.var.coord(id_reps_m_c_leg(i_reps,1),:);
    hold on;
    scatter3(foot_reps(:,1),foot_reps(:,2),foot_reps(:,3),50,'filled','MarkerFaceColor',[1 0 0]);
    hold on;
    scatter3(m_reps(:,1),m_reps(:,2),m_reps(:,3),50,'filled','MarkerFaceColor',[0 1 0]);
end
hold on;
f_test=f.int_comp.(int_name){2};
im=id_reps_m_c_leg(:,1);
quiver3(mod.mod{i_mod(2)}.var.coord(im,1),mod.mod{i_mod(2)}.var.coord(im,2),mod.mod{i_mod(2)}.var.coord(im,3),...
       f_test(im,1),f_test(im,2),f_test(im,3),'linewidth',2);
end
end

function V=potential(r1,r2,k)
   V=0.5*k*(sum((r1-r2).^2,2));
end


