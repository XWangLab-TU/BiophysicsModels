function [loc_relaxed,m] = locDyn(m,f,varargin)
%--------------------------------------------------------------------------
        % locDyn performs the local dynamics relax @ModMembrane
        % input: 
        % m - a @ModMembrane object
        % f - @TypForce object indicating forces
        % output:
        % loc_relaxed - index of relax status
        % optional:
        % see variable arguments
        %   See also TimeEval
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addParameter('mex_avail', true, @islogical);
ip.addParameter('local', 0, @isnumeric); %0: global; 1: regular local; 2: extension; 3: shrinkage
ip.addParameter('edg_exo', [], @isnumeric);%for extending or merging the one edge to cross barriers
ip.addParameter('D', 10, @isnumeric);
ip.addParameter('D2', 1000, @isnumeric);
ip.addParameter('sMax', 0.1, @isnumeric);
ip.addParameter('dt', 1e-3, @isnumeric); %1e-3
ip.addParameter('nt', 10000, @isnumeric);
ip.addParameter('mu', 1000, @isnumeric);
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(m,f,varargin{:});
%----------------------------------------------------------------------------------------
edg_exo=ip.Results.edg_exo;
sMax=ip.Results.sMax;
dt=ip.Results.dt;
nt=ip.Results.nt;
Non=numel(m.var.id_on_coord);
mex_avail=ip.Results.mex_avail;
local=ip.Results.local;
D=ip.Results.D;
plot_or_not=ip.Results.plot_or_not;
mu=ip.Results.mu;
%----------------------------------------------------------------------------------------
Vpm=m.pm.Vdh;
%----------------------------------------------------------------------------------------
if mex_avail == false
%----------------------------------------------------------------------------------------
for it=1:nt
%----------------------------------------------------------------------------------------
d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
u = d./r;

i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
i = floor(r/m.pm.dr+0.5)-i_shift;
f_edg=f.int_stored.ModMembrane.fn(i).*u;
f.int_const.ModMembrane=zeros(m.var.n_coord,3);
for i_coord = 1:Non
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)-sum(f_edg(m.var.edge_all(m.var.id_on_edg,1)==m.var.id_on_coord(i_coord),:),1);
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)+sum(f_edg(m.var.edge_all(m.var.id_on_edg,2)==m.var.id_on_coord(i_coord),:),1);
end
%----------------------------------------------------------------------------------------
id_tem1 = r<Vpm.rl_min;  id_tem2 = r>Vpm.rl_max; 
% id_tem1 = r<Vpm.r_best_min;  id_tem2 = r>Vpm.r_best_max; 
n_tem1 = length(id_tem1(id_tem1));
n_tem2 = length(id_tem2(id_tem2));
%----------------------------------------------------------------------------------------
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
%----------------------------------------------------------------------------------------
elseif local==2
    n_r=numel(r);
    id_all=1:n_r;
    id_tem3=id_all(m.var.id_on_edg==edg_exo);
    if sum(m.var.id_on_edg(id_tem2)==edg_exo) == 1
        n_tem2=n_tem2-1;
        id_tem2(id_tem3)=false;
    end
    if (r(id_tem3)>Vpm.rd_min ) && (n_tem2 == 0)
        loc_relaxed = 1;
    elseif (r(id_tem3)<Vpm.rl_max ) && (n_tem2 == 0)
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
%----------------------------------------------------------------------------------------
elseif local==3
    n_r=numel(r);
    id_all=1:n_r;
    id_tem3=id_all(m.var.id_on_edg==edg_exo);
    if sum(m.var.id_on_edg(id_tem1)==edg_exo) == 1
        n_tem1=n_tem1-1;
        id_tem1(id_tem3)=false;
    end

    if (r(id_tem3)<Vpm.rs_max ) && (n_tem1 == 0)
        loc_relaxed = 1;
    elseif (r(id_tem3)>Vpm.rl_min ) && (n_tem1 == 0)
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
%----------------------------------------------------------------------------------------
end
%----------------------------------------------------------------------------------------
    f_mag = max(sqrt(sum(f.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:).^2,2)));
    dt_tem=dt;
    if f_mag*dt_tem*m.pm.mu > sMax
        dt_tem = sMax/(f_mag*m.pm.mu);
    end
    f.int_const.ModMembrane(m.var.id_on_coord,:)=f.int_const.ModMembrane(m.var.id_on_coord,:)+f.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)*dt_tem*m.pm.mu;
end
%----------------------------------------------------------------------------------------
else
    %%
    % A:edge
    edgTem=[m.var.id_on_edg;edg_exo];
    A=m.var.edge_all(edgTem,:);
    [C,~,ic]=unique(A);
    C2=(1:numel(C))';
    A=C2(ic);
    A=reshape(A,[numel(A)/2,2]);
    
    idOrg=[m.var.id_on_coord;m.var.id_bound];
    
    [idNew,~]=sort(idOrg);
    
    [~,iNtoO] = ismember(idOrg,idNew);
    id_on_coord=iNtoO(1:Non);
    
    coord=m.var.coord(idNew,:);
    
    i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
    
    rmSeed=floor(rand(1,1)*1000+0.5);
    
%     pmc=[sMax,dt,nt,i_shift,m.pm.dr,Non,Vpm.rl_min,Vpm.rl_max,local,D,m.pm.mu,Vpm.rd_min,m.pm.Vdw.rs_max,rmSeed];
    pmc=[sMax,dt,nt,i_shift,m.pm.dr,Non,Vpm.r_best_min,Vpm.r_best_max,local,D,mu,Vpm.rd_min,m.pm.Vdh.rs_max,rmSeed];%XXX
    %----------------------------------------------------------------------
    %%
    if plot_or_not==true
    d = (coord(A(:,2),:) - coord(A(:,1),:));
    r = sqrt(sum(d.^2,2));
    u = d./r;
    i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
    i = floor(r/m.pm.dr+0.5)-i_shift;
    f_edg=f.int_stored.ModMembrane.fn(i).*u;
    n_coord=size(coord,1); 
    fc=zeros(n_coord,3);
for i_coord = 1:Non
    icOn=id_on_coord(i_coord);
    fc(icOn,:) = fc(icOn,:)-sum(f_edg(A(:,1)==icOn,:),1);
    fc(icOn,:) = fc(icOn,:)-sum(f_edg(A(:,1)==icOn,:),1);
end
    hold on;
    scatter3(coord(:,1),coord(:,2),coord(:,3),'filled');hold on;
    scatter3(coord(id_on_coord,1),coord(id_on_coord,2),coord(id_on_coord,3),'filled');hold on;
    quiver3(coord(A(:,1),1),coord(A(:,1),2),coord(A(:,1),3),f_edg(:,1),f_edg(:,2),f_edg(:,3)); hold on;
    plot3(coord(A(end,:),1),coord(A(end,:),2),coord(A(end,:),3),'linewidth',2);
    xlimAdj=[xlim,ylim,zlim];
    xlimAdj=[min(xlimAdj),max(xlimAdj)];
    xlim(xlimAdj);ylim(xlimAdj);zlim(xlimAdj);
    hold on;
    for i=1:numel(r)
        plot3([coord(A(i,2),1); coord(A(i,1),1)],[coord(A(i,2),2); coord(A(i,1),2)],[coord(A(i,2),3); coord(A(i,1),3)],'linewidth',2);hold on;
    end
%     scatter3(coord(A(end,1),1),coord(A(end,1),2),coord(A(end,1),3),'filled');hold on;
%     scatter3(coord(A(end,2),1),coord(A(end,2),2),coord(A(end,2),3),'filled');hold on;
%     quiver3(coord(A(end,1),1),coord(A(end,1),2),coord(A(end,1),3),f_edg(end,1),f_edg(end,2),f_edg(end,3));
fprintf('~~~~~ %d\n',local);
pause;
    end
    %%
    %----------------------------------------------------------------------
%     coord=coord_save;
    coord_save=coord;
    [coord,loc_relaxed,fc,check]=ModMembrane.locDynMex...
       (coord,...
        id_on_coord,...
        A,...
        pmc,...
        f.int_stored.ModMembrane.fn/m.pm.Vdh.V0*20,...
        f.int_stored.ModMembrane.R,...
        f.int_stored.ModMembrane.in,...
        f.int_stored.ModMembrane.rg);
    d = (coord(A(:,2),:) - coord(A(:,1),:));
    r = sqrt(sum(d.^2,2));
    disp([min(r),max(r)]);
    check=floor(check+0.5);
    if ((check==1) && (m.pm.remeshScheme==0))
        disp('out of boarder!');
        disp('out of boarder!');
        pause;
    end
    %%
%     scatter3(coord(:,1),coord(:,2),coord(:,3),'filled');hold on;
%     quiver3(coord(:,1),coord(:,2),coord(:,3),fc(:,1),fc(:,2),fc(:,3));
    loc_relaxed=floor(loc_relaxed+0.5);
    m.var.coord(idNew,:)=coord;
end


