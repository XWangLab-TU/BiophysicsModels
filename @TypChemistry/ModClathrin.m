function [mod,changed] = ModClathrin(ch,mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ch', @(x) isa(x,'TypChemistry'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('dt', 0.01, @isnumeric);
ip.addParameter('dr_adp', 1, @isnumeric);
ip.addParameter('lim_xyz', [-5 5;-5 5;0 5], @isnumeric);
ip.parse(ch,mod,varargin{:});
%----------------------------------------------------------------------------------------
i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane,mod.i_mod.ModMemAdapter];  %1: clathrin, 2: membrane, 3: membrane adapter
lim_xyz=ip.Results.lim_xyz;
dt=ip.Results.dt;
dr_adp=ip.Results.dr_adp;
mod=ip.Results.mod;
mex=ip.Results.mex;
if isempty(mex)
    mex_avail=false;
else
    mex_avail=true;
end
%----------------------------------------------------------------------------------------
%%
[k_on_list,n_list]=getKonList(mod,i_mod,dt);
% changed=false;
%%
if n_list>0
%%
changed=false;
%[id_mesh_mem] = ComMath.getMeshID(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_on_coord,:),mod.Mesh.coord,mod.Mesh.range,mod.Mesh.d);
id_waist=[8 16 24];
id_feet=[49 50 51];

%shift_xyz=ComMath.grid3((-1:d_xyz:1)',(-1:d_xyz:1)',(-1:d_xyz:1)');
% if mod.mod{i_mod(1)}.var.n_coord==1
%     shift_xyz=rand(100,3)*norm(mod.mod{1}.var.coord_org(1,1:3)-mod.mod{1}.var.coord_org(2,1:3))*5;
% else
     shift_xyz=rand(6,3)*norm(mod.mod{1}.var.coord_org(1,1:3)-mod.mod{1}.var.coord_org(2,1:3))*5;
% end
%shift_xyz=zeros(1,3);

n_xyz=size(shift_xyz,1);
d_all=inf(mod.mod{i_mod(1)}.pm.n_test_all,n_xyz);
i_leg_all=zeros(mod.mod{i_mod(1)}.pm.n_test_all,n_xyz);
n_on_mem=zeros(mod.mod{i_mod(1)}.pm.n_test_all,n_xyz);
i_on_mem=false(mod.mod{i_mod(1)}.pm.n_test_all,3,n_xyz);
i_c_mem=zeros(mod.mod{i_mod(1)}.pm.n_test_all,3,n_xyz);
i_leg_rand=randsample(3,3);
%%
for i_list=1:n_list
    i_c=k_on_list(i_list,1);
    i_leg_host=k_on_list(i_list,2);
%--------------------------------------------------------------------------
junction=[i_c i_leg_host];
        [~,n_ring] = getTopology(mod.mod{i_mod(1)},junction);
        if (n_ring(1) < 6) && (n_ring(2) < 6)
%%
%--------------------------------------------------------------------------
waist_host=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_waist(i_leg_host),1:3),1,mod.mod{i_mod(1)}.var.O(:,:,i_c))+mod.mod{i_mod(1)}.var.coord(i_c,:);
CtrlPt_host=mod.mod{i_mod(1)}.get_r_from_a([mod.mod{i_mod(1)}.var.ctrl_pt_org(i_leg_host,:);mod.mod{i_mod(1)}.var.ctrl_pt_org(i_leg_host+3,:)],...
                                           2,mod.mod{i_mod(1)}.var.O(:,:,i_c))+mod.mod{i_mod(1)}.var.coord(i_c,:);
%--------------------------------------------------------------------------
    for i_xyz=1:n_xyz
%--------------------------------------------------------------------------
    coord_new=waist_host+shift_xyz(i_xyz,:);
    d_check=coord_new'-lim_xyz;
    if (numel(d_check((d_check(:,1)>0),1))==3)&&(numel(d_check((d_check(:,2)<0),1))==3)
%--------------------------------------------------------------------------
%%
coord_host_hub=mod.mod{i_mod(1)}.var.coord_org(end,1:3)+mod.mod{i_mod(1)}.var.coord(i_c,:);
%r_tem1=mod.mod{i_mod(1)}.get_r_from_a([0,0,1],1,mod.mod{i_mod(1)}.var.O(:,:,i_c));
for i_test=1:mod.mod{i_mod(1)}.pm.n_test_all
   %-----------------------------------------------------------------------
   foot=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet,1:3),3,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
   %-----------------------------------------------------------------------
   waist=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_waist,1:3),3,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
   %----------------------------------------------------------------------- 
   %[id_mesh_foot] = ComMath.getMeshID(foot,mod.Mesh.coord,mod.Mesh.range,mod.Mesh.d);
   %----------------------------------------------------------------------- 
   CtrlPt=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.ctrl_pt_org,...
                                           6,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
   %-----------------------------------------------------------------------
   occupied=mod.mod{i_mod(3)}.var.occupied;
   d_all(i_test,i_xyz)=0;
   n_adp_bound=0;
   for i_foot=1:3
       d=sqrt(sum((foot(i_foot,:)-mod.mod{i_mod(3)}.var.coord).^2,2));
       d(occupied)=dr_adp+1;
       [d_min,id_min]=min(d);
       if d_min<dr_adp
       n_on_mem(i_test,i_xyz)=n_on_mem(i_test,i_xyz)+1;
       i_c_mem(i_test,i_foot,i_xyz)=id_min;
       i_on_mem(i_test,i_foot,i_xyz)=true;
       d_all(i_test,i_xyz)=d_all(i_test,i_xyz)+d_min;
       occupied(id_min)=true;
       n_adp_bound=n_adp_bound+1;
       end
   end
   if n_adp_bound<mod.mod{i_mod(1)}.pm.n_min_adp
       d_all(i_test,i_xyz)=inf;
   elseif n_adp_bound>mod.mod{i_mod(1)}.pm.n_max_adp
       d_all(i_test,i_xyz)=inf;
   else
       d_all(i_test,i_xyz)=d_all(i_test,i_xyz)/n_adp_bound;
   end
   %-----------------------------------------------------------------------
%    [id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepSingle(mod.TypForce,mod,coord_new,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test));
%    if isempty(id_reps_m_c_leg)
%        allow_add=true;
%    else
% %        id_tem = sum(abs(id_reps_m_c_leg(:,2)-i_c),2)==0;
% %        if (isempty(n_in_mem(id_tem))) || (sum(n_in_mem(id_tem))<3)
% %            allow_add=true;
% %        else
% %            allow_add=false;
% %        end
%        allow_add=false;
%    end
   allow_add=true;
   
   if  (allow_add==true) %(n_on_mem(i_test,i_xyz)>=1)
       [d_tem,i_leg_all(i_test,i_xyz)]=min(sqrt(sum((waist-coord_host_hub).^2,2))...
                                          +sqrt(sum((coord_new-waist_host).^2,2))...
                                          +sqrt(sum((CtrlPt_host(1,:)-CtrlPt(4:6,:)).^2,2))...
                                          +sqrt(sum((CtrlPt_host(2,:)-CtrlPt(1:3,:)).^2,2))); 
       d_all(i_test,i_xyz)=d_all(i_test,i_xyz)+d_tem;                                                 
   else
       d_all(i_test,i_xyz)=inf;
   end
   

%     r_tem2=mod.mod{i_mod(1)}.get_r_from_a([0,0,1],1,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test));
%     cos_tem=dot(r_tem1,r_tem2)/norm(r_tem1)/norm(r_tem2);
%     if cos_tem<mod.mod{i_mod(1)}.pm.cos_min
%         d_all(i_test-i_c)=inf;
%     end
end
%--------------------------------------------------------------------------
    end
%--------------------------------------------------------------------------
    end
%--------------------------------------------------------------------------    
%%
if numel(d_all(~isinf(d_all)))>0
[~,id_xyz]=min(d_all,[],1);
id_keep=zeros(2,1);
min_tem=inf;
for i_xyz=1:n_xyz
    if min_tem>d_all(id_xyz(i_xyz),i_xyz)
       min_tem=d_all(id_xyz(i_xyz),i_xyz);
       id_keep=[id_xyz(i_xyz),i_xyz];
    end
end
%--------------------------------------------------------------------------
coord_new=waist_host+shift_xyz(id_keep(2),:);
mod.mod{i_mod(1)}.var.coord=[mod.mod{i_mod(1)}.var.coord;coord_new];
mod.mod{i_mod(1)}.var.n_coord=mod.mod{i_mod(1)}.var.n_coord+1;
mod.mod{i_mod(1)}.var.a=[mod.mod{i_mod(1)}.var.a;mod.mod{i_mod(1)}.pm.a_test_all(id_keep(1),:)];
mod.mod{i_mod(1)}.var.ang_a.Phi=[mod.mod{i_mod(1)}.var.ang_a.Phi;mod.mod{i_mod(1)}.pm.ang_a_test_all(id_keep(1),1)];
mod.mod{i_mod(1)}.var.O=mod.mod{i_mod(1)}.Omega(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
mod.mod{i_mod(1)}.var.E=mod.mod{i_mod(1)}.getE(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
mod.mod{i_mod(1)}.var.connect(i_leg_host,1,i_c)=mod.mod{i_mod(1)}.var.n_coord;
mod.mod{i_mod(1)}.var.connect(i_leg_host,2,i_c)=i_leg_all(id_keep(1),id_keep(2));
mod.mod{i_mod(1)}.var.connect=cat(3,mod.mod{i_mod(1)}.var.connect,zeros(3,2));
mod.mod{i_mod(1)}.var.connect(i_leg_all(id_keep(1),id_keep(2)),1,mod.mod{i_mod(1)}.var.n_coord)=i_c;
mod.mod{i_mod(1)}.var.connect(i_leg_all(id_keep(1),id_keep(2)),2,mod.mod{i_mod(1)}.var.n_coord)=i_leg_host;
mod.mod{i_mod(1)}.var.bound(mod.mod{i_mod(1)}.var.n_coord,:)=i_on_mem(id_keep(1),:,id_keep(2));
%mod.mod{i_mod(1)}.var.bound(mod.mod{i_mod(1)}.var.n_coord,:)=true(1,3);
for i_leg_tem=1:3
    if i_c_mem(id_keep(1),i_leg_tem,id_keep(2)) > 0
      mod.mod{i_mod(3)}.var.id_ModClathrin(i_c_mem(id_keep(1),i_leg_tem,id_keep(2)),:)=[mod.mod{i_mod(1)}.var.n_coord i_leg_tem];
      mod.mod{i_mod(3)}.var.occupied(i_c_mem(id_keep(1),i_leg_tem,id_keep(2)))=true;
    end
end
changed=true;
end
%--------------------------------------------------------------------------
        end
%--------------------------------------------------------------------------  
   if changed==true
       break;
   end
%--------------------------------------------------------------------------
end

% fig=figure;
% plot(mod.mod{i_mod},'col',[],'f',fig);
% hold on;
% r_tem1=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c));
% quiver3(mod.mod{i_mod}.var.coord(i_c,1),mod.mod{i_mod}.var.coord(i_c,2),mod.mod{i_mod}.var.coord(i_c,3),r_tem1(1),r_tem1(2),r_tem1(3),...
%         10,'linewidth',2);
% hold on;
% r_tem2=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c+1));
% quiver3(mod.mod{i_mod}.var.coord(i_c+1,1),mod.mod{i_mod}.var.coord(i_c+1,2),mod.mod{i_mod}.var.coord(i_c+1,3),r_tem2(1),r_tem2(2),r_tem2(3),...
%         10,'linewidth',2);
end
end
function [k_on_list,n_list]=getKonList(mod,i_mod,dt)
    k_on_list=[];
for i_c=1:mod.mod{i_mod(1)}.var.n_coord
    i_leg_host_all=randsample(3,3);
    for i_leg_host_try=1:3
        i_leg_host=i_leg_host_all(i_leg_host_try);
        if (rand(1,1) < dt*mod.mod{i_mod(1)}.pm.k_on) && (mod.mod{i_mod(1)}.var.connect(i_leg_host,1,i_c) == 0) 
            k_on_list=[k_on_list;[i_c i_leg_host]];
            break;
        end
    end
end
n_list=size(k_on_list,1);
id_tem=randsample(n_list,n_list);
k_on_list=k_on_list(id_tem,:);
end