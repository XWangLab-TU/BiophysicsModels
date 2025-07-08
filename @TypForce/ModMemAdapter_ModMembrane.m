function [f,loc_relaxed,abnormal_length,m] = ModMemAdapter_ModMembrane(f,m,adp,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('adp', @(x) isa(x,'ModMemAdapter'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('f_const_only', false, @islogical);
ip.addParameter('f_comp_only', false, @islogical);
ip.addParameter('local', 0, @isnumeric); %0: global; 1: regular local; 2: extension; 3: shrinkage
ip.addParameter('edg_exo', 0, @isnumeric);%for extending or merging the one edge to cross barriers
ip.addParameter('D', 10000, @isnumeric);
ip.addParameter('D2', 1000, @isnumeric);
ip.addParameter('Pflat', false, @islogical);
ip.addParameter('add_spring', false, @islogical); %additional force at unsplitable/unmergable vertices
ip.addParameter('split', [], @isstruct);
ip.addParameter('merge', [], @isstruct);
ip.parse(f,m,adp,varargin{:});
%----------------------------------------------------------------------------------------
mex=ip.Results.mex;
if isempty(mex)
    mex_avail=false;
else
    mex_avail=true;
end
f_const_only=ip.Results.f_const_only;
f_comp_only=ip.Results.f_comp_only;
local=ip.Results.local;
D=ip.Results.D;
D2=ip.Results.D2;
edg_exo=ip.Results.edg_exo;
split=ip.Results.split;
merge=ip.Results.merge;
%----------------------------------------------------------------------------------------
if isempty(f.int_stored.ModMembrane)
    f.int_stored.ModMembrane = ModMembrane8Const(f,m);
end
loc_relaxed=false;
abnormal_length=false;
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
if f_comp_only==false
d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
u = d./r;
i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
i = floor(r/m.pm.dr+0.5)-i_shift;
f_edg=f.int_stored.ModMembrane.fn(i).*u;
f.int_const.ModMembrane=zeros(m.var.n_coord,3);

for i_coord = 1:numel(m.var.id_on_coord)
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)-sum(f_edg(m.var.edge_all(m.var.id_on_edg,1)==m.var.id_on_coord(i_coord),:),1);
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)+sum(f_edg(m.var.edge_all(m.var.id_on_edg,2)==m.var.id_on_coord(i_coord),:),1);
end
id_tem1 = r<m.pm.Vdw.rl_min;  id_tem2 = r>m.pm.Vdw.rl_max; 
n_tem1 = length(id_tem1(id_tem1));
n_tem2 = length(id_tem2(id_tem2));
if (n_tem1 == 0) && (n_tem2 == 0)
    abnormal_length=false;
else
    abnormal_length=true;
end
if local==1
    if (n_tem1 == 0) && (n_tem2 == 0)
    loc_relaxed = true;
    else
    if n_tem1 > 0
        f_rand = randsample(f.int_stored.ModMembrane.R,n_tem1).*u(id_tem1,:)*D;%+randn(n_tem1,3)*D2;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:)-f_rand;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:)+f_rand;
    end
    if n_tem2 > 0
        f_rand = -randsample(f.int_stored.ModMembrane.R,n_tem2).*u(id_tem2,:)*D;%+randn(n_tem2,3)*D2;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:)-f_rand;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:)+f_rand;
    end
    end
elseif local==2
    n_r=numel(r);
    id_all=1:n_r;
    id_tem3=id_all(m.var.id_on_edg==edg_exo);
    if sum(m.var.id_on_edg(id_tem2)==edg_exo) == 1
        n_tem2=n_tem2-1;
        id_tem2(id_tem3)=false;
    end

    if (r(id_tem3)>m.pm.Vdw.rd_min ) && (n_tem2 == 0)
        loc_relaxed = 1;
    elseif (r(id_tem3)<m.pm.Vdw.rl_max ) && (n_tem2 == 0)
        loc_relaxed = 2;
    else
        if n_tem1 > 0
            f_rand = randsample(f.int_stored.ModMembrane.R,n_tem1).*u(id_tem1,:)*D;%+randn(n_tem1,3)*D2;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:)-f_rand;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:)+f_rand;
        end
        if n_tem2 > 0
            f_rand = -randsample(f.int_stored.ModMembrane.R,n_tem2).*u(id_tem2,:)*D;%+randn(n_tem2,3)*D2;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:)-f_rand;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:)+f_rand;
        end
        f_rand = randsample(f.int_stored.ModMembrane.R,1).*u(id_tem3,:)*D;%+randn(n_tem1,3)*D2;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),1),:)-f_rand;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),2),:)+f_rand;
    end
elseif local==3
    n_r=numel(r);
    id_all=1:n_r;
    id_tem3=id_all(m.var.id_on_edg==edg_exo);
    if sum(m.var.id_on_edg(id_tem1)==edg_exo) == 1
        n_tem1=n_tem1-1;
        id_tem1(id_tem3)=false;
    end
    %%
%     fig=figure('units','normalized','outerposition',[0 0 0.5 1]);
%     plot(m,'linestyle','-','f',fig,'FaceAlpha',1,'light_pos',[0.5 0.6 0.8]);
%     hold on;
%     id_edg_test=m.var.id_on_edg;
%     n_test=numel(id_edg_test);
%     for i=1:n_test
%     plot3([m.var.coord(m.var.edge_all(id_edg_test(i),1),1),m.var.coord(m.var.edge_all(id_edg_test(i),2),1)],...
%       [m.var.coord(m.var.edge_all(id_edg_test(i),1),2),m.var.coord(m.var.edge_all(id_edg_test(i),2),2)],...
%       [m.var.coord(m.var.edge_all(id_edg_test(i),1),3),m.var.coord(m.var.edge_all(id_edg_test(i),2),3)],'linewidth',3,'color',[1 1 0]);hold on;
%     end
%     id_edg_test=edg_exo;
%     plot3([m.var.coord(m.var.edge_all(id_edg_test,1),1),m.var.coord(m.var.edge_all(id_edg_test,2),1)],...
%       [m.var.coord(m.var.edge_all(id_edg_test,1),2),m.var.coord(m.var.edge_all(id_edg_test,2),2)],...
%       [m.var.coord(m.var.edge_all(id_edg_test,1),3),m.var.coord(m.var.edge_all(id_edg_test,2),3)],'linewidth',5,'color',[1 0 0]);hold on;
%     zlim([-10 10]);
    %%
    if (r(id_tem3)<m.pm.Vdw.rs_max ) && (n_tem1 == 0)
        loc_relaxed = 1;
    elseif (r(id_tem3)>m.pm.Vdw.rl_min ) && (n_tem1 == 0)
        loc_relaxed = 2;
    else
        if n_tem1 > 0
            f_rand = randsample(f.int_stored.ModMembrane.R,n_tem1).*u(id_tem1,:)*D;%+randn(n_tem1,3)*D2;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),1),:)-f_rand;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem1),2),:)+f_rand;
        end
        if n_tem2 > 0
            f_rand = -randsample(f.int_stored.ModMembrane.R,n_tem2).*u(id_tem2,:)*D;%+randn(n_tem2,3)*D2;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),1),:)-f_rand;
            f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem2),2),:)+f_rand;
        end
        f_rand = -randsample(f.int_stored.ModMembrane.R,1).*u(id_tem3,:)*D;%+randn(n_tem1,3)*D2;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),1),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),1),:)-f_rand;
        f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),2),:) = f.int_const.ModMembrane(m.var.edge_all(m.var.id_on_edg(id_tem3),2),:)+f_rand;
    end
end
end
%----------------------------------------------------------------------------------------
if f_const_only==false
    if m.pm.mu==0
        f.int_rand.ModMemAdapter_ModMembrane = zeros(m.var.n_coord,3);
    else
        f.int_rand.ModMemAdapter_ModMembrane = randn(m.var.n_coord,3)*sqrt(2*m.pm.mu*m.pm.kBT*m.pm.dt)/m.pm.dt/m.pm.mu;
    end
if mex_avail == false
%%
[m] = f.ModMemAdapter_ModMembrane8Helfrich(m,adp,'idx',m.var.id_on_coord);
[m] = getUface(m);
[m] = ModMembrane8uK(f,m);
Cadp = getC(adp,m);
m = f.ModMemAdapter_ModMembrane8Helfrich(m,adp,'idx',m.var.id_on_coord,'Cadp',Cadp,'init',false);
coord_org = m.var.coord;
H_org = m.var.f.H;
K_org=m.var.f.K;
kH_org = m.var.f.kH;
A_org = m.var.f.A;
f.int_comp.ModMemAdapter_ModMembrane = zeros(m.var.n_coord,3);
m.var.f.dH = m.var.f.dH*0;
dr_H=0.00001;
for i_on = 1:m.var.n_on_coord
    i=m.var.id_on_coord(i_on);
    for k=1:3
       m.var.coord = coord_org; m.var.coord(i,k) = m.var.coord(i,k) + dr_H;
       idx = [i,m.var.j_T(i,1:m.var.val(i))];
       [~,id_tem]=intersect(idx,m.var.id_on_coord);
       idx=idx(id_tem);
       m = f.ModMemAdapter_ModMembrane8Helfrich(m,adp,'C',Cadp,'idx',idx,'init',false);
       H_new = m.var.f.H;
%        id_tem=1:m.var.n_coord;
%        for i_tem=1:max(size(idx))
%            id_tem(id_tem==idx(i_tem))=[];
%        end
%        norm(kH_org(id_tem,:) - m.var.f.kH(id_tem,:))
%        pause
       [min_tem,id_min_tem]=min(abs(i-adp.var.id_ModMembrane));
       if min_tem==0
           if adp.var.kind(id_min_tem)==1
               f.int_comp.ModMemAdapter_ModMembrane(i,k) = -m.pm.k_c*(H_new-H_org)/dr_H;
           elseif adp.var.kind(id_min_tem)==2
               f.int_comp.ModMemAdapter_ModMembrane(i,k) = -adp.pm.k_c*(H_new-H_org)/dr_H;
           else
               error('wrong kind');
           end
       else
           f.int_comp.ModMemAdapter_ModMembrane(i,k) = -m.pm.k_c*(H_new-H_org)/dr_H;
       end
       m.var.f.dH(i,k) = (H_new-H_org)/dr_H;
       m.var.f.kH(idx)=kH_org(idx);
       m.var.f.A(idx)=A_org(idx);
       m.var.f.K(idx,:)=K_org(idx,:);
    end
end


m.var.coord = coord_org;
%% test plot
% fig=figure('units','normalized','outerposition',[0 0 1 1]); 
% view_theta=0; view_phi=50; x_lim=[-15 15];
% 
% plot(m,'linestyle','-','f',fig,'FaceAlpha',1,'view_ang',[view_theta view_phi]);
% plot(adp,m,'f',fig,'plot_clathrin_adp',false);
% xlim(x_lim);ylim(x_lim);zlim([-3 3]);
% 
% id_tem=Cadp>0.1;
% scatter3(m.var.coord(id_tem,1),m.var.coord(id_tem,2),m.var.coord(id_tem,3),40,'filled','markerfacecolor',[1 0 0]);
% id_tem=Cadp<-0.1;
% scatter3(m.var.coord(id_tem,1),m.var.coord(id_tem,2),m.var.coord(id_tem,3),40,'filled','markerfacecolor',[0 1 0]);
% 
% %f_test=f.int_comp.ModMemAdapter_ModMembrane;
% %f_test=m.var.f.u_K;
% K_tem = m.var.f.K./m.var.f.A; f_test = K_tem./sqrt(sum(K_tem.*K_tem,2));
% hold on;
% quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),f_test(:,1),f_test(:,2),f_test(:,3),'linewidth',2);
%--------------------------------------------------------------------------
else
    pmc=zeros(7,1);
    pmc(1) = m.pm.nt;
    pmc(2) = m.pm.dt;
    pmc(3) = m.pm.P;
    pmc(4) = m.pm.k_c;
    pmc(5) = m.pm.k_e;
    pmc(6) = m.pm.dr;
    pmc(7) = adp.pm.k_c;
    
    j_T = m.var.j_T; j_T(isnan(j_T)) = 0;
    
    Cadp=zeros(m.var.n_coord,1);
    Cadp(adp.var.id_ModMembrane)=adp.var.C;
    
    Kadp=m.pm.k_c*ones(m.var.n_coord,1);
    Kadp(adp.var.id_ModMembrane(adp.var.kind==2))=adp.pm.k_c;

    [f_b,f_p,u_K,kH,f_AV]=mex.ModMemAdapter_ModMembrane...
       (m.var.coord,...
        pmc,...
        m.var.edge_all,...
        m.var.face_unq,...
        j_T,...
        m.var.id_on_coord,...
        m.var.n_node',...
        m.var.T_s',...
        m.var.T_e',...
        Cadp,...
        Kadp);
%% test plot
% f_test=f_b;
% hold on;
% quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),f_test(:,1),f_test(:,2),f_test(:,3),'linewidth',2);
%%
    m.var.f.kH=kH;
    m.var.f.u_K=u_K;
    m.var.f.f_AV=f_AV;
    if ip.Results.Pflat==true
    f_p(m.var.coord(:,3)<0.0001,:)=0;
    end
    f.int_comp.ModMemAdapter_ModMembrane=f_b+f_p+f_AV;
    id_edg_all=(1:m.var.n_edg)';
    if ip.Results.add_spring==true
        i_no_split=id_edg_all(~split.can);
        n_no_split=numel(i_no_split);
        i_no_merge=id_edg_all(~merge.can);
        n_no_merge=numel(i_no_merge);
        d = (m.var.coord(m.var.edge_all(i_no_split,2),:) - m.var.coord(m.var.edge_all(i_no_split,1),:));
        r = sqrt(sum(d.^2,2));
        u = d./r;
        f_add=(m.pm.k_const*(1./tanh((m.pm.Vdw.rl_max-r)*m.pm.Vdw.k_w))-1).*u;
        for i=1:n_no_split
            f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_split(i),1),:)=f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_split(i),1),:)+f_add(i,:);
            f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_split(i),2),:)=f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_split(i),2),:)-f_add(i,:);
        end
        d = (m.var.coord(m.var.edge_all(i_no_merge,2),:) - m.var.coord(m.var.edge_all(i_no_merge,1),:));
        r = sqrt(sum(d.^2,2));
        u = d./r;
        f_add=(m.pm.k_const*(1./tanh((r-m.pm.Vdw.rl_min)*m.pm.Vdw.k_w))-1).*u;
        for i=1:n_no_merge
            f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_merge(i),1),:)=f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_merge(i),1),:)-f_add(i,:);
            f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_merge(i),2),:)=f.int_comp.ModMemAdapter_ModMembrane(m.var.edge_all(i_no_merge(i),2),:)+f_add(i,:);
        end
    end
end
end

