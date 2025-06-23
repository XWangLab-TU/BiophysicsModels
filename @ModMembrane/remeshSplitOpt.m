function [m,remeshed,edg_add] = remeshSplitOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim,varargin)
%--------------------------------------------------------------------------
        % remeshSplitOpt performs the splitting on @ModMembrane,
        % input: 
        % m - @ModMembrane input
        % j - flipping index
        % k - flipping index
        % id_ring_edg - ring id near splitting edge
        % edg_add_org - merging edge
        % i_edg - index of original edge
        % rLim - acceptable range of edge length
        % output:
        % remeshed - true: remesh done, false: remesh gave up
        % edg_add - id of changed edges
        % optional:
        % see variable arguments
        %   See also remeshFlipOpt,remeshMergeOpt
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
var_new.coord = [m.var.coord;0.5*(m.var.coord(m.var.edge_all(i_edg,1),:)+m.var.coord(m.var.edge_all(i_edg,2),:))];
var_new.n_coord = var_new.n_coord+1;
%----------------------------------------------------------------------
i_n = m.var.n_coord+1; i_e = m.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
dens_new=var_new.dens;
densTem=sum(dens_new(i_e))/3;
dens_new(i_e)=dens_new(i_e)/3*2;
dens_new=[dens_new; densTem];
%----------------------------------------------------------------------
edg_rem = [i_edg;id_ring_edg(sum(m.var.edge_all(id_ring_edg,1)==j(1),2)+sum(m.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
    id_ring_edg(sum(m.var.edge_all(id_ring_edg,1)==j(2),2)+sum(m.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
id_rem = false(m.var.n_edg,1);
id_rem(edg_rem) = true;
var_new.edge_all(id_rem,:) = [];
r_tem = [sqrt(sum((var_new.coord(edg_add_org(:,1),:)-var_new.coord(edg_add_org(:,2),:)).^2,2)),...
    sqrt(sum((var_new.coord(edg_add_org(:,3),:)-var_new.coord(edg_add_org(:,4),:)).^2,2)),...
    sqrt(sum((var_new.coord(edg_add_org(:,5),:)-var_new.coord(edg_add_org(:,6),:)).^2,2)),...
    sqrt(sum((var_new.coord(edg_add_org(:,7),:)-var_new.coord(edg_add_org(:,8),:)).^2,2))];
[~,id_sort] = sort(std(r_tem-m.pm.l0,[],2));
n_try=numel(id_sort);
edg_save=edg_add_org;
var_new_save=var_new;
successTem=false;
mSave=m;
for i_try=1:n_try
    id_tem=id_sort(i_try);
    edg_tem = edg_save(id_tem,:);
    edg_tem = [edg_tem(1:2);edg_tem(3:4);edg_tem(5:6);edg_tem(7:8)];
    edg_add = [edg_tem;[m.var.edge_all(i_edg,1) i_n];[i_n m.var.edge_all(i_edg,2)];[i_n j(1)];[i_n j(2)]];
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
    face_add = zeros(8,3);
    i_f = 1;
    if sum(edg_add(i_f,:) == i_e(1),2)+sum(edg_add(i_f,:) == i_e(2),2) == 1
        face_add(i_f,:) = [i_e(1) i_n j(1)]; face_add(i_f+1,:) = [j(1) k(1) i_e(1)];
    else
        face_add(i_f,:) = [i_e(1) i_n k(1)]; face_add(i_f+1,:) = [i_n j(1) k(1)];
    end
    i_f = 2;
    if sum(edg_add(i_f,:) == i_e(1),2)+sum(edg_add(i_f,:) == i_e(2),2) == 1
        face_add(i_f*2-1,:) = [i_n i_e(2) j(1)]; face_add(i_f*2,:) = [i_e(2) k(2) j(1)];
    else
        face_add(i_f*2-1,:) = [i_n i_e(2) k(2)]; face_add(i_f*2,:) = [k(2) j(1) i_n];
    end
    i_f = 3;
    if sum(edg_add(i_f,:) == i_e(1),2)+sum(edg_add(i_f,:) == i_e(2),2) == 1
        face_add(i_f*2-1,:) = [i_e(2) i_n j(2)]; face_add(i_f*2,:) = [i_e(2) j(2) k(3)];
    else
        face_add(i_f*2-1,:) = [i_e(2) i_n k(3)]; face_add(i_f*2,:) = [k(3) i_n j(2)];
    end
    i_f = 4;
    if sum(edg_add(i_f,:) == i_e(1),2)+sum(edg_add(i_f,:) == i_e(2),2) == 1
        face_add(i_f*2-1,:) = [i_n i_e(1) j(2)]; face_add(i_f*2,:) = [j(2) i_e(1) k(4)];
    else
        face_add(i_f*2-1,:) = [i_n i_e(1) k(4)]; face_add(i_f*2,:) = [i_n k(4) j(2)];
    end
    var_new.face_unq = [var_new.face_unq;face_add];
    %--------------------------------------------------------------------------
    [m.var,topologicalDefect] = m.remeshAddVertex(m.pm,m.var,var_new.coord,var_new.edge_all,var_new.face_unq,...
        'dens',dens_new);
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
