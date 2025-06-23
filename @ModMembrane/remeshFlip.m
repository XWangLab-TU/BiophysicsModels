function [flip] = remeshFlip(m,M,idTooLong,varargin)
%--------------------------------------------------------------------------
        % flip performs the flipping manipulation on @ModMembrane,
        % preparing data needed for remesh
        % input: 
        % m - @ModMembrane input
        % M - @model input
        % idTooLong - index of edges too long that need flipping
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
ip.addRequired('idTooLong', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('only_V', false, @islogical);
ip.parse(m,M,idTooLong,varargin{:});
%--------------------------------------------------------------------------------------------------------
plot_or_not=ip.Results.plot_or_not;
%--------------------------------------------------------------------------
            i_mod=M.i_mod.ModMembrane;  %membrane
%--------------------------------------------------------------------------------------------------------
%%
var=struct('n_edg',M.mod{i_mod}.var.n_edg,...
           'face_unq',M.mod{i_mod}.var.face_unq,...
           'edge_all',M.mod{i_mod}.var.edge_all,...
           'val',M.mod{i_mod}.var.val,...
           'n_ver',M.mod{i_mod}.var.n_coord);

flip = struct('can',[],'edgFlip',[],'j',[],'id_ring_edg',[]);
flip.edgFlip = cell(var.n_edg,1); flip.j = cell(var.n_edg,1); flip.id_ring_edg = cell(var.n_edg,1);
flip.can = true(var.n_edg,1);
%%
%==========================================================================
n_flip=numel(idTooLong);
id_ring_edg_save=cell(var.n_edg,1);
for i_sm = 1:n_flip
    i_edg=idTooLong(i_sm);

if isempty(id_ring_edg_save{i_edg})
    [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
    id_ring_edg_save{i_edg}=id_ring_edg;
else
    id_ring_edg=id_ring_edg_save{i_edg};
end
if numel(intersect(id_ring_edg,M.mod{i_mod}.var.id_on_edg))==numel(id_ring_edg) 
%--------------------------------------------------------------------------    
flip.id_ring_edg{i_edg} = id_ring_edg;
id_tem = sum(var.face_unq == var.edge_all(i_edg,1),2)+sum(var.face_unq == var.edge_all(i_edg,2),2);
id_tem = var.face_unq(id_tem==2,:);
id_tem2 = id_tem(1,:);
j = zeros(2,1);
j(1) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));
id_tem2 = id_tem(2,:);
j(2) = id_tem2(id_tem2 ~= var.edge_all(i_edg,1) & id_tem2 ~= var.edge_all(i_edg,2));

flip.j{i_edg} = j;
%==========================================================================
%after flipping, add 2 tri and remove 1 at old vertex, 
%remove 2 and add 1 at new vertex
if M.mod{i_mod}.var.val(var.edge_all(i_edg,1))-1<M.mod{i_mod}.pm.n_val_min
    flip.can(i_edg)=false;
elseif M.mod{i_mod}.var.val(var.edge_all(i_edg,2))-1<M.mod{i_mod}.pm.n_val_min
    flip.can(i_edg)=false;
elseif M.mod{i_mod}.var.val(j(1))+1>M.mod{i_mod}.pm.n_val_max
    flip.can(i_edg)=false;
elseif M.mod{i_mod}.var.val(j(2))+1>M.mod{i_mod}.pm.n_val_max
    flip.can(i_edg)=false;
end
%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------------------------------------------------
flip.edgFlip = [j(1) j(2)]; 
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------topological
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
    scatter3(coord(flip.j{i_edg}(1),1),coord(flip.j{i_edg}(1),2),coord(flip.j{i_edg}(1),3),20,'filled','markerfacecolor',[0 1 0]);
    hold on;
    scatter3(coord(flip.j{i_edg}(2),1),coord(flip.j{i_edg}(2),2),coord(flip.j{i_edg}(2),3),20,'filled','markerfacecolor',[1 1 0]);
    hold on;
    
    id_ring_edg=flip.id_ring_edg{i_edg};
    n_ring_edg=numel(id_ring_edg);
    for i=1:n_ring_edg
        hold on;
        plot3([coord(edg(id_ring_edg(i),1),1),coord(edg(id_ring_edg(i),2),1)],[coord(edg(id_ring_edg(i),1),2),coord(edg(id_ring_edg(i),2),2)],[coord(edg(id_ring_edg(i),1),3),coord(edg(id_ring_edg(i),2),3)],...
          'linewidth',2,'color',[0 0 1]);
    end
    xlim(x_lim);ylim(x_lim);zlim(x_lim);
end
