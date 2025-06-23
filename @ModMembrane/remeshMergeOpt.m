function [m,remeshed] = remeshMergeOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim,varargin)
%--------------------------------------------------------------------------
        % remeshMergeOpt performs the Merging on @ModMembrane,
        % input: 
        % m - @ModMembrane input
        % j - flipping index
        % k - flipping index
        % id_ring_edg - ring id near merging edge
        % edg_add_org - merging edge
        % i_edg - index of original edge
        % rLim - acceptable range of edge length
        % output:
        % remeshed - true: remesh done, false: remesh gave up
        % optional:
        % see variable arguments
        %   See also remeshFlipOpt,remeshSplitOpt
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('j', @(x) isnumeric(x));
ip.addRequired('k', @(x) isnumeric(x));
ip.addRequired('id_ring_edg', @(x) isnumeric(x));
ip.addRequired('edg_add_org', @(x) isnumeric(x));
ip.addRequired('i_edg', @(x) isnumeric(x));
ip.addRequired('rLim', @(x) isnumeric(x));
ip.parse(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim,varargin{:});
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%%
%----------------------------------------------------------------------
var_new=m.var;
i = m.var.edge_all(i_edg,:);
var_new.coord(i(1),:) = 0.5*(m.var.coord(i(1),:)+m.var.coord(i(2),:));
%----------------------------------------------------------------------
dens_new=var_new.dens;
dens_new(i(1))=dens_new(i(1))+dens_new(i(2));
dens_new(i(2))=[];
%----------------------------------------------------------------------
i_e = m.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
edg_rem = [i_edg;...
    id_ring_edg(sum(m.var.edge_all(id_ring_edg,1)==j(1),2)+sum(m.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
    id_ring_edg(sum(m.var.edge_all(id_ring_edg,1)==j(2),2)+sum(m.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
id_rem = false(m.var.n_edg,1);
id_rem(edg_rem) = true;
var_new.edge_all(id_rem,:) = [];
try
r_tem = [sqrt(sum((var_new.coord(edg_add_org(:,1),:)-var_new.coord(edg_add_org(:,2),:)).^2,2)),...
    sqrt(sum((var_new.coord(edg_add_org(:,3),:)-var_new.coord(edg_add_org(:,4),:)).^2,2))];
catch
    disp('wrong');
    disp('wrong');
end
[~,id_sort] = sort(std(r_tem-m.pm.l0,[],2));
n_try=numel(id_sort);
successTem=false;
edg_save=edg_add_org;
var_new_save=var_new;
mSave=m;
for i_try=1:n_try
    id_tem=id_sort(i_try);
    edg_tem = edg_save(id_tem,:);
    edg_tem = [edg_tem(1:2);edg_tem(3:4)];
    edg_add = edg_tem;
    var_new.edge_all = [var_new.edge_all;edg_add];
    %----------------------------------------------------------------------
    face_rem = zeros(6,1);
    n_face = size(m.var.face_unq,1);
    id_all = 1:n_face;
    try
        i_tem = sum(m.var.face_unq==i_e(1),2) + sum(m.var.face_unq==j(1),2) + sum(m.var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
    catch
        disp('wrong');
        disp('wrong');
    end
    i_tem = sum(m.var.face_unq==i_e(1),2) + sum(m.var.face_unq==i_e(2),2) + sum(m.var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
    i_tem = sum(m.var.face_unq==i_e(2),2) + sum(m.var.face_unq==k(2),2) + sum(m.var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
    i_tem = sum(m.var.face_unq==i_e(2),2) + sum(m.var.face_unq==j(2),2) + sum(m.var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
    i_tem = sum(m.var.face_unq==i_e(1),2) + sum(m.var.face_unq==i_e(2),2) + sum(m.var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
    i_tem = sum(m.var.face_unq==i_e(1),2) + sum(m.var.face_unq==j(2),2) + sum(m.var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
    id_rem = false(n_face,1);
    id_rem(face_rem) = true;
    var_new.face_unq(id_rem,:) = [];
    %--------------------------------------------------------------------------
    face_add = zeros(4,3);
    i_f = 1;
    if sum(edg_add(i_f,:) == i_e(1),2) == 1
        face_add(i_f,:) = [j(1) i_e(1) k(2)]; face_add(i_f+1,:) = [j(1) k(1) i_e(1)];
    else
        face_add(i_f,:) = [j(1) k(1) k(2)]; face_add(i_f+1,:) = [k(1) i_e(1) k(2)];
    end
    i_f = 2;
    if sum(edg_add(i_f,:) == i_e(1),2) == 1
        face_add(i_f*2-1,:) = [j(2) k(3) i_e(1)]; face_add(i_f*2,:) = [i_e(1) k(4) j(2)];
    else
        face_add(i_f*2-1,:) = [j(2) k(3) k(4)]; face_add(i_f*2,:) = [k(4) k(3) i_e(1)];
    end
    
    var_new.face_unq = [var_new.face_unq;face_add];
    
    var_new.face_unq(var_new.face_unq==i_e(2)) = i_e(1);
    id_tem = var_new.face_unq>i_e(2);
    var_new.face_unq(id_tem) = var_new.face_unq(id_tem)-1;
    
    var_new.edge_all(var_new.edge_all==i_e(2)) = i_e(1);
    id_tem = var_new.edge_all>i_e(2);
    var_new.edge_all(id_tem) = var_new.edge_all(id_tem)-1;
    
    edg_add(edg_add>i_e(2))=edg_add(edg_add>i_e(2))-1;
    
    var_new.coord(i(2),:) = [];
    var_new.n_coord = var_new.n_coord-1;
    
    [m.var,topologicalDefect] = m.remeshAddVertex(m.pm,m.var,var_new.coord,var_new.edge_all,var_new.face_unq,'dens',dens_new);
    if topologicalDefect==true
        m.failInfo='topologicalDefect';
        successTem=false;
    else
        r = sqrt(sum(([m.var.coord(m.var.edge_all(:,2),1),m.var.coord(m.var.edge_all(:,2),2),m.var.coord(m.var.edge_all(:,2),3)]...
            -[m.var.coord(m.var.edge_all(:,1),1),m.var.coord(m.var.edge_all(:,1),2),m.var.coord(m.var.edge_all(:,1),3)]).^2,2));
        
        if (max(r)<rLim(2)) && (min(r)>rLim(1))
            successTem=true;
            break;
        elseif (i_try==n_try)
            successTem=false;
        end
    end
    m=mSave;
    var_new=var_new_save;
end

if successTem==true
    [m] = getUface(m);
    remeshed=true;
else
    m=mSave;
    remeshed=false;
end
                
end
