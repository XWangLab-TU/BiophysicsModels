function [var, id_s, id_m] = femmConfChange(pm,var,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('var', @(x) isstruct(x));
ip.parse(pm, var,varargin{:});
%--------------------------------------------------------------------------------------------------------
changed = true; 
relaxed = true;
i_try = 0;
while (changed == true)
%====================================================================================================================
%====================================================================================================================
%%
[split,merge] = femmSplitMerge(pm,var);
r = sqrt(sum(([var.ver(var.edge_all(:,2),1),var.ver(var.edge_all(:,2),2),var.ver(var.edge_all(:,2),3)]...
             -[var.ver(var.edge_all(:,1),1),var.ver(var.edge_all(:,1),2),var.ver(var.edge_all(:,1),3)]).^2,2));

id_split = r>pm.Vdw.rl_max;
id_merge = r<pm.Vdw.rl_min;

id_all = (1:var.n_edg)';

id_split = id_all(id_split);
id_merge = id_all(id_merge);

id_sm = [id_split;id_merge];
s_or_m = [ones(length(id_split),1);zeros(length(id_merge),1)];
n_sm = size(id_sm,1);
%%
if i_try == 0
    id_s=id_split;
    id_s=id_s(split.can(id_s));
    id_m=id_merge;
    id_m=id_m(merge.can(id_m));
end
i_try = i_try+1;
%%
if n_sm > 0
i_sm = randsample(n_sm,1);
i_edg = id_sm(i_sm);
%fprintf('%d, %d, %d\n',i_edg, s_or_m(i_sm), var.n_ver);
if s_or_m(i_sm) == 1
%====================================================================================================================   
if split.can(i_edg)
    var_new=var;
    var_new.ver = [var.ver;0.5*(var.ver(var.edge_all(i_edg,1),:)+var.ver(var.edge_all(i_edg,2),:))];
    var_new.n_ver = var_new.n_ver+1;
    %----------------------------------------------------------------------
    j = split.j{i_edg}; k = split.k{i_edg}; id_ring_edg = split.id_ring_edg{i_edg}; i_n = var.n_ver+1; i_e = var.edge_all(i_edg,:);
    %----------------------------------------------------------------------
    edg_rem = [i_edg;id_ring_edg(sum(var.edge_all(id_ring_edg,1)==j(1),2)+sum(var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
                     id_ring_edg(sum(var.edge_all(id_ring_edg,1)==j(2),2)+sum(var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
    id_rem = false(var.n_edg,1);
    id_rem(edg_rem) = true;
    var_new.edge_all(id_rem,:) = [];
    edg_tem = split.edg_add{i_edg};
    r_tem = [sqrt(sum((var_new.ver(edg_tem(:,1),:)-var_new.ver(edg_tem(:,2),:)).^2,2)),...
             sqrt(sum((var_new.ver(edg_tem(:,3),:)-var_new.ver(edg_tem(:,4),:)).^2,2)),...
             sqrt(sum((var_new.ver(edg_tem(:,5),:)-var_new.ver(edg_tem(:,6),:)).^2,2)),...
             sqrt(sum((var_new.ver(edg_tem(:,7),:)-var_new.ver(edg_tem(:,8),:)).^2,2))];
    [~,id_tem] = min(std(r_tem-pm.l0,[],2));
    edg_tem = edg_tem(id_tem,:);
    edg_tem = [edg_tem(1:2);edg_tem(3:4);edg_tem(5:6);edg_tem(7:8)];
    edg_add = [edg_tem;[var.edge_all(i_edg,1) i_n];[i_n var.edge_all(i_edg,2)];[i_n j(1)];[i_n j(2)]];
    var_new.edge_all = [var_new.edge_all;edg_add];
    %----------------------------------------------------------------------
    face_rem = zeros(6,1);
    n_face = size(var.face_unq,1);
    id_all = 1:n_face;
    i_tem = sum(var.face_unq==i_e(1),2) + sum(var.face_unq==j(1),2) + sum(var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
    i_tem = sum(var.face_unq==i_e(1),2) + sum(var.face_unq==i_e(2),2) + sum(var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
    i_tem = sum(var.face_unq==i_e(2),2) + sum(var.face_unq==k(2),2) + sum(var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
    i_tem = sum(var.face_unq==i_e(2),2) + sum(var.face_unq==j(2),2) + sum(var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
    i_tem = sum(var.face_unq==i_e(1),2) + sum(var.face_unq==i_e(2),2) + sum(var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
    i_tem = sum(var.face_unq==i_e(1),2) + sum(var.face_unq==j(2),2) + sum(var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
    id_rem = false(n_face,1);
    id_rem(face_rem) = true;
    var_new.face_unq(id_rem,:) = [];
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    [var] = femmAddVertex(pm,var,var_new.ver,var_new.edge_all,var_new.face_unq);
    %plotMem(var.ver,var.face_unq,gcf,ones(var.n_ver,3),'FaceAlpha', 1, 'LineStyle','-');
    [var,relaxed] = femmLocRelax(pm,var,edg_add);
%     if relaxed == false
%         fprintf('local relax failed\n');
%         %changed = false;
%     end
    %----------------------------------------------------------------------
else
    %----------------------------------------------------------------------
    [var,relaxed] = femmLocRelax(pm,var,var.edge_all(i_edg,:));
%     if relaxed == false
%         fprintf('local relax failed\n');
%         %changed = false;
%     end
    %----------------------------------------------------------------------
end
%====================================================================================================================
else
%====================================================================================================================    
if merge.can(i_edg)
    %----------------------------------------------------------------------
    edge_all_new = var.edge_all; id_all = 1:max(size(edge_all_new));
    face_new = var.face_unq;
    i = var.edge_all(i_edg,:);
    i = sort(i);
    id_tem = sum(var.face_unq == i(1),2) + sum(var.face_unq == i(2),2);
    id_tem = id_tem == 2;
    j = var.face_unq(id_tem,:);
    j((j==i(1))|(j==i(2))) = [];
    face_new(id_tem,:) = [];
    id_tem = face_new>i(2);
    face_new(face_new==i(2)) = i(1);
    face_new(id_tem) = face_new(id_tem)-1;   
    id_tem1 = sum(var.edge_all == i(1),2) + sum(var.edge_all == i(2),2);
    id_tem1 = id_tem1==2;
    id_tem2 = sum(var.edge_all == i(2),2) + sum(var.edge_all == j(1),2);
    id_tem2 = id_tem2==2;
    id_tem3 = sum(var.edge_all == i(2),2) + sum(var.edge_all == j(2),2);
    id_tem3 = id_tem3==2;
    id_tem4 = var.edge_all == i(2);
    edge_all_new(id_tem4) = i(1); 
    edge_all_new(id_tem1 | id_tem2 | id_tem3,:) = [];
    id_tem = edge_all_new>i(2);
    edge_all_new(id_tem) = edge_all_new(id_tem)-1;
    ver_new = var.ver;
    ver_new(i(2),:) = [];
    ver_new(i(1),:) = 0.5*(var.ver(i(1),:)+var.ver(i(2),:));
    [var] = femmAddVertex(pm,var,ver_new,edge_all_new,face_new);
    %----------------------------------------------------------------------
%     plotMem(var.ver,var.face_unq,gcf,ones(var.n_ver,3),'FaceAlpha', 1, 'LineStyle','-'); hold on;
%     scatter3(var.ver(i(1),1),var.ver(i(1),2),var.ver(i(1),3),40,'filled'); 
    %----------------------------------------------------------------------
    edg_add = var.edge_all(sum(var.edge_all == i(1),2)>0,:);
    [var,relaxed] = femmLocRelax(pm,var,edg_add);
%     if relaxed == false
%         fprintf('local relax failed\n');
%         %changed = false;
%     end
else
    [var,relaxed] = femmLocRelax(pm,var,var.edge_all(i_edg,:));
%     if relaxed == false
%         fprintf('local relax failed\n');
%         %changed = false;
%     end
end
%====================================================================================================================
end
else
    changed = false;
end
%====================================================================================================================
%====================================================================================================================
end

 