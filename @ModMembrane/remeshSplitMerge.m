function [split,merge] = remeshSplitMerge(m,M,id_split,id_merge,varargin)
%--------------------------------------------------------------------------
        % flip performs the splitting/merging manipulation on @ModMembrane,
        % preparing data needed for remesh
        % input: 
        % m - @ModMembrane input
        % M - @model input
        % id_split - index of edges too long that need splitting
        % id_merge - index of edges too short that need merging
        % output:
        % flip - data structure of whether an edge can be flipped and etc.
        % optional:
        % see variable arguments
        %   See also remesh, remeshSplitMerge
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('id_split', @(x) isnumeric(x));
ip.addRequired('id_merge', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('only_V', false, @islogical);
ip.addParameter('ignoreCanOrNot', false, @islogical);
ip.parse(m,M,id_split,id_merge,varargin{:});
%--------------------------------------------------------------------------------------------------------
plot_or_not=ip.Results.plot_or_not;
ignoreCanOrNot=ip.Results.ignoreCanOrNot;
%--------------------------------------------------------------------------
            i_mod=M.i_mod.ModMembrane;  %membrane
%--------------------------------------------------------------------------------------------------------
%%
var=struct('n_edg',M.mod{i_mod}.var.n_edg,...
           'face_unq',M.mod{i_mod}.var.face_unq,...
           'edge_all',M.mod{i_mod}.var.edge_all,...
           'val',M.mod{i_mod}.var.val,...
           'n_ver',M.mod{i_mod}.var.n_coord);

split = struct('can',[],'edg_add',[],'j',[],'k',[],'id_ring_edg',[]);
split.edg_add = cell(var.n_edg,1); split.j = cell(var.n_edg,1); split.k = cell(var.n_edg,1); split.id_ring_edg = cell(var.n_edg,1);
split.can = true(var.n_edg,1);
merge = struct('can',[],'edg_add',[],'j',[],'k',[],'id_ring_edg',[]);
merge.edg_add = cell(var.n_edg,1); merge.j = cell(var.n_edg,1); merge.k = cell(var.n_edg,1); merge.id_ring_edg = cell(var.n_edg,1);
merge.can = true(var.n_edg,1);
%%
%==========================================================================
n_split=numel(id_split);
n_merge=numel(id_merge);
id_ring_edg_save=cell(var.n_edg,1);
for i_sm = 1:n_split+n_merge
    if i_sm<=n_split
        i_edg=id_split(i_sm);
    else
        i_edg=id_merge(i_sm-n_split);
    end

if isempty(id_ring_edg_save{i_edg})
    [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
    id_ring_edg_save{i_edg}=id_ring_edg;
else
    id_ring_edg=id_ring_edg_save{i_edg};
end
if numel(intersect(id_ring_edg,M.mod{i_mod}.var.id_on_edg))==numel(id_ring_edg) 
%--------------------------------------------------------------------------    
split.id_ring_edg{i_edg} = id_ring_edg;
merge.id_ring_edg{i_edg} = id_ring_edg;
id_tem = sum(var.face_unq == var.edge_all(i_edg,1),2)+sum(var.face_unq == var.edge_all(i_edg,2),2);
id_tem = var.face_unq(id_tem==2,:);
id_tem2 = id_tem(1,:);
[~,i_tem1] = min(abs(id_tem2-var.edge_all(i_edg,1)));
[~,i_tem2] = min(abs(id_tem2-var.edge_all(i_edg,2)));
i_tem1 = i_tem1 + 1;
if i_tem1 > 3
    i_tem1 = 1;
end
j = zeros(2,1);
if i_tem1 == i_tem2
    j(1) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));
    id_tem2 = id_tem(2,:);
    j(2) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));
else
    id_tem2 = id_tem(2,:);
    j(1) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));
    id_tem2 = id_tem(1,:);
    j(2) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));
end
split.j{i_edg} = j;
merge.j{i_edg} = j;
%--------------------------------------------------------------------------
k = zeros(4,1);
id_tem = sum(var.face_unq == var.edge_all(i_edg,1),2)+sum(var.face_unq == j(1),2);
id_tem = var.face_unq(id_tem==2,:);
k(1) = id_tem(id_tem~=var.edge_all(i_edg,1) & id_tem~=var.edge_all(i_edg,2) & id_tem ~= j(1));
id_tem = sum(var.face_unq == var.edge_all(i_edg,2),2)+sum(var.face_unq == j(1),2);
id_tem = var.face_unq(id_tem==2,:);
k(2) = id_tem(id_tem~=var.edge_all(i_edg,1) & id_tem~=var.edge_all(i_edg,2) & id_tem ~= j(1));
id_tem = sum(var.face_unq == var.edge_all(i_edg,2),2)+sum(var.face_unq == j(2),2);
id_tem = var.face_unq(id_tem==2,:);
k(3) = id_tem(id_tem~=var.edge_all(i_edg,1) & id_tem~=var.edge_all(i_edg,2) & id_tem ~= j(2));
id_tem = sum(var.face_unq == var.edge_all(i_edg,1),2)+sum(var.face_unq == j(2),2);
id_tem = var.face_unq(id_tem==2,:);
k(4) = id_tem(id_tem~=var.edge_all(i_edg,1) & id_tem~=var.edge_all(i_edg,2) & id_tem ~= j(2));
split.k{i_edg} = k;
merge.k{i_edg} = k;
%==========================================================================
%--------------------------------------------------------------------------
i_n = var.n_ver+1;
edg_add = ones(16,8); val_tem = [var.val(var.edge_all(i_edg,1:2))'-2,4,var.val(j)'-1,var.val(k)'];
val_new = val_tem.*ones(16,1);
%--------------------------------------------------------------------------
edg_add(1:8,1) = k(1); edg_add(1:8,2) = i_n;
edg_add(9:16,1) = var.edge_all(i_edg,1); edg_add(9:16,2) = j(1);
val_new(1:8,6) = val_new(1:8,6)+1; val_new(1:8,3) = val_new(1:8,3)+1;
val_new(9:16,1) = val_new(9:16,1)+1; val_new(9:16,4) = val_new(9:16,4)+1;
%--------------------------------------------------------------------------
edg_add([1:4,9:12]',3) = k(2); edg_add([1:4,9:12]',4) = i_n;
edg_add([5:8,13:16]',3) = var.edge_all(i_edg,2); edg_add([5:8,13:16]',4) = j(1);
val_new([1:4,9:12]',7) = val_new([1:4,9:12]',7)+1; val_new([1:4,9:12]',3) = val_new([1:4,9:12]',3)+1;
val_new([5:8,13:16]',2) = val_new([5:8,13:16]',2)+1; val_new([5:8,13:16]',4) = val_new([5:8,13:16]',4)+1;
%--------------------------------------------------------------------------
edg_add([1:2,5:6,9:10,13:14]',5) = k(3); edg_add([1:2,5:6,9:10,13:14]',6) = i_n;
edg_add([3:4,7:8,11:12,15:16]',5) = var.edge_all(i_edg,2); edg_add([3:4,7:8,11:12,15:16]',6) = j(2);
val_new([1:2,5:6,9:10,13:14]',8) = val_new([1:2,5:6,9:10,13:14]',8)+1; val_new([1:2,5:6,9:10,13:14]',3) = val_new([1:2,5:6,9:10,13:14]',3)+1;
val_new([3:4,7:8,11:12,15:16]',2) = val_new([3:4,7:8,11:12,15:16]',2)+1; val_new([3:4,7:8,11:12,15:16]',5) = val_new([3:4,7:8,11:12,15:16]',5)+1;
%--------------------------------------------------------------------------
edg_add(1:2:16,7) = k(4); edg_add(1:2:16,8) = i_n;
edg_add(2:2:16,7) = var.edge_all(i_edg,1); edg_add(2:2:16,8) = j(2);
val_new(1:2:16,9) = val_new(1:2:16,9) +1; val_new(1:2:16,3) = val_new(1:2:16,3) +1;
val_new(2:2:16,1) = val_new(2:2:16,1)+1; val_new(2:2:16,5) = val_new(2:2:16,5)+1;
%--------------------------------------------------------------------------
can_split = true(16,1);
can_split(sum((val_new>M.mod{i_mod}.pm.n_val_max) | (val_new<M.mod{i_mod}.pm.n_val_min),2)>0) = false;
if isempty(can_split(can_split==true))
    split.can(i_edg) = false;
else
    split.edg_add{i_edg} = edg_add(can_split==true,:);
end
if ignoreCanOrNot==true
    split.can(i_edg) = true;
    split.edg_add{i_edg} = edg_add;
end
%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------------------------------------------------
edg_add = ones(4,4); 
val_tem = [sum(var.val(var.edge_all(i_edg,1:2)))-8+2,var.val(j)'-2,var.val(k)'];
val_new = val_tem.*ones(4,1);
%--------------------------------------------------------------------------
edg_add(1:2,1)=k(1); edg_add(1:2,2)=k(2);
edg_add(3:4,1)=var.edge_all(i_edg,1); edg_add(3:4,2) = j(1);
val_new(1:2,4) = val_new(1:2,4)+1; val_new(1:2,5) = val_new(1:2,5)+1;
val_new(3:4,1) = val_new(3:4,1)+1; val_new(3:4,2) = val_new(3:4,2)+1;
%--------------------------------------------------------------------------
edg_add(1,3)=k(3); edg_add(1,4)=k(4);
val_new(1,6) = val_new(1,6)+1; val_new(1,7) = val_new(1,7)+1;
edg_add(2,3)=var.edge_all(i_edg,1); edg_add(2,4)=j(2);
val_new(2,1) = val_new(2,1)+1; val_new(2,3) = val_new(2,3)+1;
%--------------------------------------------------------------------------
edg_add(3,3)=k(3); edg_add(3,4)=k(4);
val_new(3,6) = val_new(3,6)+1; val_new(3,7) = val_new(3,7)+1;
edg_add(4,3)=var.edge_all(i_edg,1); edg_add(4,4)=j(2);
val_new(4,1) = val_new(4,1)+1; val_new(4,3) = val_new(4,3)+1;
%--------------------------------------------------------------------------
can_merge = true(4,1);
can_merge(sum((val_new>M.mod{i_mod}.pm.n_val_max) | (val_new<M.mod{i_mod}.pm.n_val_min),2)>0) = false;
if isempty(can_merge(can_merge==true))
    merge.can(i_edg) = false;
else
    merge.edg_add{i_edg} = edg_add(can_merge==true,:);
end
if ignoreCanOrNot==true
    merge.can(i_edg) = true;
    merge.edg_add{i_edg} = edg_add;
end
%--------------------------------------------------------------------------
%==========================================================================
else
    split.can(i_edg) = false;
    merge.can(i_edg) = false;
end
%--------------------------------------------------------------------------topological
%constraint: forbid merging from closing an open triangular neck (failed)
% if merge.can(i_edg)==true
%     i_e=mod.mod{i_mod}.var.edge_all(i_edg,:);
%     id_tem=sum(mod.mod{i_mod}.var.edge_all==i_e(1),2)+sum(mod.mod{i_mod}.var.edge_all==i_e(2),2);
%     id_tem=id_tem==1;
%     id_nb=mod.mod{i_mod}.var.edge_all(id_tem,:);
%     id_nb=unique(id_nb);
%     id_nb(id_nb==i_e(1))=[];
%     id_nb(id_nb==i_e(2))=[];
%     n_nb=numel(id_nb);
%     for i_nb=1:n_nb
%         id_tem=sum(mod.mod{i_mod}.var.edge_all==i_e(1),2)...
%                +sum(mod.mod{i_mod}.var.edge_all==i_e(2),2)...
%               +sum(mod.mod{i_mod}.var.edge_all==id_nb(i_nb),2);
%         edge_tem=mod.mod{i_mod}.var.edge_all(id_tem==2,:);
%         if size(edge_tem,1)>2
%             merge.can(i_edg)=false;
%             break;
%         end
%     end
% % patch('Vertices',m.var.coord,'Faces',m.var.face_unq(id_tem,:),'FaceVertexCData',rand(m.var.n_coord,3),'FaceColor','interp','facealpha',1);
% % hold on;
% % plot3([m.var.coord(i_e(1),1),m.var.coord(i_e(2),1)],...
% % [m.var.coord(i_e(1),2),m.var.coord(i_e(2),2)],...
% % [m.var.coord(i_e(1),3),m.var.coord(i_e(2),3)],'linewidth',3);
% %--------------------------------------------------------------------------
% end
%--------------------------------------------------------------------------
end
%==========================================================================

%==========================================================================
%%
if plot_or_not==true
    fig=figure('units','normalized','outerposition',[0 0 0.5 1]); 
    x_lim=[-15 15];
    plot(M.mod{i_mod},'linestyle','-','f',fig,'FaceAlpha',1);
    hold on;
    coord=M.mod{i_mod}.var.coord;
    edg=M.mod{i_mod}.var.edge_all;
    plot3([coord(edg(i_edg,1),1),coord(edg(i_edg,2),1)],[coord(edg(i_edg,1),2),coord(edg(i_edg,2),2)],[coord(edg(i_edg,1),3),coord(edg(i_edg,2),3)],...
          'linewidth',2,'color',[1 0 0]);
    hold on;
    scatter3(coord(merge.j{i_edg}(1),1),coord(merge.j{i_edg}(1),2),coord(merge.j{i_edg}(1),3),20,'filled','markerfacecolor',[0 1 0]);
    hold on;
    scatter3(coord(merge.j{i_edg}(2),1),coord(merge.j{i_edg}(2),2),coord(merge.j{i_edg}(2),3),20,'filled','markerfacecolor',[1 1 0]);
    hold on;
    scatter3(coord(merge.k{i_edg}(1),1),coord(merge.k{i_edg}(1),2),coord(merge.k{i_edg}(1),3),20,'filled','markerfacecolor',[1 0 0]);
    hold on;
    scatter3(coord(merge.k{i_edg}(2),1),coord(merge.k{i_edg}(2),2),coord(merge.k{i_edg}(2),3),20,'filled','markerfacecolor',[1 0 1]);
    hold on;
    scatter3(coord(merge.k{i_edg}(3),1),coord(merge.k{i_edg}(3),2),coord(merge.k{i_edg}(3),3),20,'filled','markerfacecolor',[1 0 1]);
    hold on;
    scatter3(coord(merge.k{i_edg}(4),1),coord(merge.k{i_edg}(4),2),coord(merge.k{i_edg}(4),3),20,'filled','markerfacecolor',[1 0 1]);
    n_ring_edg=numel(id_ring_edg);
    for i=1:n_ring_edg
        hold on;
        plot3([coord(edg(id_ring_edg(i),1),1),coord(edg(id_ring_edg(i),2),1)],[coord(edg(id_ring_edg(i),1),2),coord(edg(id_ring_edg(i),2),2)],[coord(edg(id_ring_edg(i),1),3),coord(edg(id_ring_edg(i),2),3)],...
          'linewidth',2,'color',[0 0 1]);
    end
    xlim(x_lim);ylim(x_lim);zlim(x_lim);
end
