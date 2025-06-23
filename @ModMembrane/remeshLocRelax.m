function [M,loc_relaxed] = remeshLocRelax(m,M,edg_add, varargin)
%--------------------------------------------------------------------------
        % flip performs the local relaxing after remeshing on @ModMembrane,
        % including the edge indicated and its 1-ring neighbors
        % input: 
        % m - @ModMembrane input
        % M - @model input
        % edg_add - index of edge to be relaxed around
        % output:
        % loc_relaxed - various cases indicating the relax results
        % optional:
        % see variable arguments
        %   See also remesh
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('edg_add', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('dt', 0.001, @isnumeric);
ip.addParameter('n_step', 100000, @isnumeric);
ip.addParameter('D', 100, @isnumeric);
ip.addParameter('local', 0, @isnumeric); %0: global; 1: regular local; 2: extension; 3: shrinkage
ip.addParameter('ring_ord', 1, @isnumeric);
ip.parse(m,M, edg_add, varargin{:});
%--------------------------------------------------------------------------------------------------------
plot_or_not = ip.Results.plot_or_not;
local=ip.Results.local;
ring_ord=ip.Results.ring_ord;
%--------------------------------------------------------------------------
            i_mod=M.i_mod.ModMembrane;  %membrane
%--------------------------------------------------------------------------
%%
edg_add_org=edg_add;
n_add = size(edg_add,1);
id_ring_ver = [];
id_ring_edg = [];
id_all = 1:M.mod{i_mod}.var.n_edg;
id_edg_add = zeros(n_add,1);
for i = 1:n_add
    i_edg = sum(M.mod{i_mod}.var.edge_all==edg_add(i,1),2)+sum(M.mod{i_mod}.var.edge_all==edg_add(i,2),2);
    i_edg = id_all(i_edg==2);
    id_edg_add(i) = i_edg;
    if ring_ord==1
        [id_ring_ver_tem,id_ring_edg_tem] = remeshRing(M.mod{i_mod},i_edg,'ring_ord',ring_ord);
    elseif ring_ord==2
        [~,~,id_ring_ver_tem,id_ring_edg_tem] = remeshRing(M.mod{i_mod},i_edg,'ring_ord',ring_ord);
    else
        error('ring_ord is wrong');
    end
    id_ring_ver = [id_ring_ver; id_ring_ver_tem];
    id_ring_edg = [id_ring_edg; id_ring_edg_tem];
    %%---------------------------------------------------------------------
%     coord=M.mod{i_mod}.var.coord;
%     edg=M.mod{i_mod}.var.edge_all;
%     scatter3(coord(id_ring_ver,1),coord(id_ring_ver,2),coord(id_ring_ver,3),'filled');hold on;
%     scatter3(coord(edg_add,1),coord(edg_add,2),coord(edg_add,3),'filled');hold on;
%     for iii=1:numel(id_ring_edg)
%     plot3(coord(edg(id_ring_edg(iii),:),1),coord(edg(id_ring_edg(iii),:),2),coord(edg(id_ring_edg(iii),:),3));hold on;
%     end
    %%---------------------------------------------------------------------
end
for i = 1:n_add
    id_ring_edg(id_ring_edg==id_edg_add(i)) = [];
    id_ring_ver(id_ring_ver==edg_add(i,1) | id_ring_ver==edg_add(i,2)) = [];
end
id_ring_edg = unique(id_ring_edg);
id_ring_edg = [id_ring_edg; id_edg_add];
id_ring_ver = unique(id_ring_ver);
n_ring_edg = size(id_ring_edg,1);
n_ring_ver = size(id_ring_ver,1);
id_all = 1:M.mod{i_mod}.var.n_coord;
id_ring_in=unique(id_all(M.mod{i_mod}.var.edge_all(id_ring_edg,:)));
for i=1:n_ring_ver
    id_ring_in(id_ring_in==id_ring_ver(i))=[];
end

if local>1
    edg_add=[edg_add;M.mod{i_mod}.var.edge_all(id_ring_edg,:)];
    edg_add = unique(edg_add,'row');
    n_add = size(edg_add,1);
id_ring_ver = [];
id_ring_edg = [];
id_all = 1:M.mod{i_mod}.var.n_edg;
id_edg_add = zeros(n_add,1);
for i = 1:n_add
    i_edg = sum(M.mod{i_mod}.var.edge_all==edg_add(i,1),2)+sum(M.mod{i_mod}.var.edge_all==edg_add(i,2),2);
    i_edg = id_all(i_edg==2);
    id_edg_add(i) = i_edg;
    if ring_ord==1
        [id_ring_ver_tem,id_ring_edg_tem] = remeshRing(M.mod{i_mod},i_edg,'ring_ord',ring_ord);
    elseif ring_ord==2
        [~,~,id_ring_ver_tem,id_ring_edg_tem] = remeshRing(M.mod{i_mod},i_edg,'ring_ord',ring_ord);
    else
        error('ring_ord is wrong');
    end
    id_ring_ver = [id_ring_ver; id_ring_ver_tem];
    id_ring_edg = [id_ring_edg; id_ring_edg_tem];
end
for i = 1:n_add
    id_ring_edg(id_ring_edg==id_edg_add(i)) = [];
    id_ring_ver(id_ring_ver==edg_add(i,1) | id_ring_ver==edg_add(i,2)) = [];
end
id_ring_edg = unique(id_ring_edg);
id_ring_edg = [id_ring_edg; id_edg_add];
id_ring_ver = unique(id_ring_ver);
n_ring_edg = size(id_ring_edg,1);
n_ring_ver = size(id_ring_ver,1);
id_all = 1:M.mod{i_mod}.var.n_coord;
id_ring_in=unique(id_all(M.mod{i_mod}.var.edge_all(id_ring_edg,:)));
for i=1:n_ring_ver
    id_ring_in(id_ring_in==id_ring_ver(i))=[];
end
end

id_on_edg_save=M.mod{i_mod}.var.id_on_edg;
n_on_edg_save=M.mod{i_mod}.var.n_on_edg;
id_on_coord_save=M.mod{i_mod}.var.id_on_coord;
n_on_coord_save=M.mod{i_mod}.var.n_on_coord;
id_bound_save=M.mod{i_mod}.var.id_bound;
n_bound_save=M.mod{i_mod}.var.n_bound;
id_all = 1:M.mod{i_mod}.var.n_edg;
if local >1
%    edg_add_org_flip=[edg_add_org(2) edg_add_org(1)];
    id_tem=sum(abs(M.mod{i_mod}.var.edge_all-edg_add_org),2)==0;
    edg_exo=id_all(id_tem);
else
    edg_exo=[];
end
M.mod{i_mod}.var.id_on_edg=id_ring_edg;
M.mod{i_mod}.var.id_on_coord=id_ring_in;
M.mod{i_mod}.var.id_bound=id_ring_ver;
[id_int,i_int]=intersect(M.mod{i_mod}.var.id_on_coord,id_bound_save);
if ~isempty(id_int)
    M.mod{i_mod}.var.id_on_coord(i_int)=[];
    M.mod{i_mod}.var.id_bound=[M.mod{i_mod}.var.id_bound;id_int];
end
%%
[loc_relaxed,M.mod{i_mod}] = locDyn(M.mod{i_mod},M.TypForce,'edg_exo',edg_exo,'nt',ip.Results.n_step,'local',local,'D',M.mod{i_mod}.pm.DlocRelax);

M.mod{i_mod}.var.id_on_edg=id_on_edg_save;
M.mod{i_mod}.var.n_on_edg=n_on_edg_save;
M.mod{i_mod}.var.id_on_coord=id_on_coord_save;
M.mod{i_mod}.var.n_on_coord=n_on_coord_save;
M.mod{i_mod}.var.id_bound=id_bound_save;
M.mod{i_mod}.var.n_bound=n_bound_save;

%%
if (plot_or_not == true) %|| (relaxed == false)
    %save('temp.mat','var','edg_add');
    m=M.mod{i_mod};
% m.var.id_on_edg=id_ring_edg;
% m.var.id_on_coord=id_ring_in;
% m.var.id_bound=id_ring_ver;

plot(m,'FaceAlpha', 1, 'LineStyle','none'); xlabel('x');ylabel('y');zlabel('z');hold on;
scatter3(m.var.coord(edg_add(:,1),1),m.var.coord(edg_add(:,1),2),m.var.coord(edg_add(:,1),3),40,'filled','MarkerFaceColor',[0 1 0]); hold on;
scatter3(m.var.coord(edg_add(:,2),1),m.var.coord(edg_add(:,2),2),m.var.coord(edg_add(:,2),3),40,'filled','MarkerFaceColor',[0 1 0]); hold on;
scatter3(m.var.coord(m.var.id_bound,1),m.var.coord(m.var.id_bound,2),m.var.coord(m.var.id_bound,3),40,'filled','MarkerFaceColor',[0 1 1]); hold on;
scatter3(m.var.coord(m.var.id_on_coord,1),m.var.coord(m.var.id_on_coord,2),m.var.coord(m.var.id_on_coord,3),40,'filled','MarkerFaceColor',[0 0.2 1]); hold on;

if local>1
    plot3([m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,1),1),m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,2),1)],...
      [m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,1),2),m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,2),2)],...
      [m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,1),3),m.var.coord(M.mod{i_mod}.var.edge_all(edg_exo,2),3)],'linewidth',5,'color',[1 0 0]);hold on;
end

for i = 1:n_add
    plot3([m.var.coord(edg_add(i,1),1),m.var.coord(edg_add(i,2),1)],...
      [m.var.coord(edg_add(i,1),2),m.var.coord(edg_add(i,2),2)],...
      [m.var.coord(edg_add(i,1),3),m.var.coord(edg_add(i,2),3)],'linewidth',3,'color',[0 1 0]);hold on;
end

for i = 1:n_ring_edg
    plot3([m.var.coord(m.var.edge_all(m.var.id_on_edg(i),1),1),m.var.coord(m.var.edge_all(m.var.id_on_edg(i),2),1)],...
      [m.var.coord(m.var.edge_all(m.var.id_on_edg(i),1),2),m.var.coord(m.var.edge_all(m.var.id_on_edg(i),2),2)],...
      [m.var.coord(m.var.edge_all(m.var.id_on_edg(i),1),3),m.var.coord(m.var.edge_all(m.var.id_on_edg(i),2),3)],'linewidth',2,'color',[1 1 0]);hold on;
end
end
