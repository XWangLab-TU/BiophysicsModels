function [mod,id_s,id_m,remeshed] = ModMemAdapter_ModMembrane(ch,mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ch', @(x) isa(x,'TypChemistry'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('f_const_only', true, @islogical);
ip.addParameter('print_or_not', false, @islogical);
ip.parse(ch,mod,varargin{:});
%----------------------------------------------------------------------------------------
f_const_only=ip.Results.f_const_only;
print_or_not=ip.Results.print_or_not;
%--------------------------------------------------------------------------
%%
            i_mod=[mod.i_mod.ModMembrane mod.i_mod.ModMemAdapter];  %1. membrane 2. ModMemAdapter
%----------------------------------------------------------------------------------------
%%
changed = true; 
remeshed=false;
i_try = 0;
while (changed == true)
    %%
    r = sqrt(sum(([mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,2),1),mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,2),2),mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,2),3)]...
        -[mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,1),1),mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,1),2),mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(:,1),3)]).^2,2));
    
    id_split = r>mod.mod{i_mod(1)}.pm.Vdw.rl_max;
    id_merge = r<mod.mod{i_mod(1)}.pm.Vdw.rl_min;
    
    id_all = (1:mod.mod{i_mod(1)}.var.n_edg)';
    
    id_split = id_all(id_split);
    id_merge = id_all(id_merge);
    [split,merge] = ModMembrane8SplitMerge(ch,mod,id_split,id_merge);

    id_sm = [id_split;id_merge];
    s_or_m = [ones(length(id_split),1);zeros(length(id_merge),1)];
    n_sm = size(id_sm,1);
    if print_or_not==true
    fprintf('abnormal edge: %d, long: %d, short: %d\n',n_sm,numel(id_split),numel(id_merge));
    end
    for i=1:n_sm
        i_e = mod.mod{i_mod(1)}.var.edge_all(id_sm(i),:);
        id_tem1=mod.mod{i_mod(2)}.var.id_ModMembrane==i_e(1);
        id_tem2=mod.mod{i_mod(2)}.var.id_ModMembrane==i_e(2);
        is_occupied1=(mod.mod{i_mod(2)}.var.occupied(id_tem1)==true);
        if isempty(is_occupied1)
            is_occupied1=false;
        elseif numel(is_occupied1)>1
            is_occupied1=sum(is_occupied1==true);
        end
        is_occupied2=(mod.mod{i_mod(2)}.var.occupied(id_tem2)==true);
        if isempty(is_occupied2)
            is_occupied2=false;
        elseif numel(is_occupied2)>1
            is_occupied2=sum(is_occupied2==true);
        end
        if (is_occupied1==true) && (is_occupied2==true)
            merge.can(id_sm(i))=false;
        end
    end
    %%
%----------------------------------------------------------------------------------------    
    if i_try == 0
    id_s=id_split;
    id_s=id_s(split.can(id_s));
    id_m=id_merge;
    id_m=id_m(merge.can(id_m));
    end
    i_try = i_try+1;
%----------------------------------------------------------------------------------------
    id_sm_can=[id_split(split.can(id_split));id_merge(merge.can(id_merge))];
    n_sm_can=numel(id_sm_can);
    id_sm_rest=[id_split(~split.can(id_split));id_merge(~merge.can(id_merge))];
    n_sm_rest=n_sm-n_sm_can;
%----------------------------------------------------------------------------------------    
    if n_sm > 0
        if n_sm_can>0
            i_sm = randsample(n_sm_can,1);
            i_edg = id_sm_can(i_sm);
        else
            i_sm = randsample(n_sm_rest,1);
            i_edg = id_sm_rest(i_sm);
        end
% i_sm = randsample(n_sm,1);
% i_edg = id_sm(i_sm);
        if s_or_m(i_sm) == 1
%----------------------------------------------------------------------------------------
            if split.can(i_edg)
%----------------------------------------------------------------------
               [mod,relaxed] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',2);
               if relaxed==2
                [mod,~] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',1);
                if print_or_not==true
                    fprintf('failed: initial extension\n');
                end
               elseif relaxed==1
%----------------------------------------------------------------------
                remeshed=true;
                var_new=mod.mod{i_mod(1)}.var;
                var_new.coord = [mod.mod{i_mod(1)}.var.coord;0.5*(mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(i_edg,1),:)+mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.edge_all(i_edg,2),:))];
                var_new.n_coord = var_new.n_coord+1;
%----------------------------------------------------------------------
                j = split.j{i_edg}; k = split.k{i_edg}; id_ring_edg = split.id_ring_edg{i_edg}; i_n = mod.mod{i_mod(1)}.var.n_coord+1; i_e = mod.mod{i_mod(1)}.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
                edg_rem = [i_edg;id_ring_edg(sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,1)==j(1),2)+sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
                    id_ring_edg(sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,1)==j(2),2)+sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
                id_rem = false(mod.mod{i_mod(1)}.var.n_edg,1);
                id_rem(edg_rem) = true;
                var_new.edge_all(id_rem,:) = [];
                edg_tem = split.edg_add{i_edg};
                r_tem = [sqrt(sum((var_new.coord(edg_tem(:,1),:)-var_new.coord(edg_tem(:,2),:)).^2,2)),...
                    sqrt(sum((var_new.coord(edg_tem(:,3),:)-var_new.coord(edg_tem(:,4),:)).^2,2)),...
                    sqrt(sum((var_new.coord(edg_tem(:,5),:)-var_new.coord(edg_tem(:,6),:)).^2,2)),...
                    sqrt(sum((var_new.coord(edg_tem(:,7),:)-var_new.coord(edg_tem(:,8),:)).^2,2))];
                [~,id_tem] = min(std(r_tem-mod.mod{i_mod(1)}.pm.l0,[],2));
                edg_tem = edg_tem(id_tem,:);
                edg_tem = [edg_tem(1:2);edg_tem(3:4);edg_tem(5:6);edg_tem(7:8)];
                edg_add = [edg_tem;[mod.mod{i_mod(1)}.var.edge_all(i_edg,1) i_n];[i_n mod.mod{i_mod(1)}.var.edge_all(i_edg,2)];[i_n j(1)];[i_n j(2)]];
                var_new.edge_all = [var_new.edge_all;edg_add];
%----------------------------------------------------------------------
                face_rem = zeros(6,1);
                n_face = size(mod.mod{i_mod(1)}.var.face_unq,1);
                id_all = 1:n_face;
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
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
                [mod.mod{i_mod(1)}.var] = ch.ModMembrane8AddVertex(mod.mod{i_mod(1)}.pm,mod.mod{i_mod(1)}.var,var_new.coord,var_new.edge_all,var_new.face_unq);
                [mod.mod{i_mod(1)}] = getUface(mod.mod{i_mod(1)});
                [mod,relaxed] = ModMembrane8LocRelax(ch,mod,edg_add,'f_const_only',f_const_only,'local',1);
               else
                   if print_or_not==true
                   fprintf('failed: initial extension reached loop limit\n');
                   end
               end
            else
%--------------------------------------------------------------------------
                [mod,relaxed] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',1);
                if print_or_not==true
                fprintf('failed: unsplittable\n');
                end
%--------------------------------------------------------------------------
            end
%----------------------------------------------------------------------------------------
        else
            if merge.can(i_edg)
%----------------------------------------------------------------------------------------
               [mod,relaxed] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',3);
               if relaxed==2
                [mod,~] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',1);
                if print_or_not==true
                fprintf('failed: initial merging\n');
                end
               elseif relaxed==1
%----------------------------------------------------------------------------------------
                remeshed=true;
                var_new=mod.mod{i_mod(1)}.var;
                i = mod.mod{i_mod(1)}.var.edge_all(i_edg,:);
                var_new.coord(i(1),:) = 0.5*(mod.mod{i_mod(1)}.var.coord(i(1),:)+mod.mod{i_mod(1)}.var.coord(i(2),:));
%----------------------------------------------------------------------
                j = merge.j{i_edg}; k = merge.k{i_edg}; id_ring_edg = merge.id_ring_edg{i_edg}; i_e = mod.mod{i_mod(1)}.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
                edg_rem = [i_edg;...
                           id_ring_edg(sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,1)==j(1),2)+sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
                           id_ring_edg(sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,1)==j(2),2)+sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
%                 edg_replace = id_ring_edg(sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,1)==i_e(2),2)+sum(mod.mod{i_mod(1)}.var.edge_all(id_ring_edg,2)==i_e(2),2) > 0);
%                 [~,id_duplicated]=intersect(edg_replace,edg_rem);
%                 edg_replace(id_duplicated)=[];
%                 edg_tem=var_new.edge_all(edg_replace,:);
%                 edg_tem(edg_tem==i_e(2))=i_e(1);
%                 var_new.edge_all(edg_replace,:)=edg_tem;
                id_rem = false(mod.mod{i_mod(1)}.var.n_edg,1);
                id_rem(edg_rem) = true;
                var_new.edge_all(id_rem,:) = [];
                edg_tem = merge.edg_add{i_edg};
                r_tem = [sqrt(sum((var_new.coord(edg_tem(:,1),:)-var_new.coord(edg_tem(:,2),:)).^2,2)),...
                    sqrt(sum((var_new.coord(edg_tem(:,3),:)-var_new.coord(edg_tem(:,4),:)).^2,2))];
                [~,id_tem] = min(std(r_tem-mod.mod{i_mod(1)}.pm.l0,[],2));
                edg_tem = edg_tem(id_tem,:);
                edg_tem = [edg_tem(1:2);edg_tem(3:4)];
                edg_add = edg_tem;
                var_new.edge_all = [var_new.edge_all;edg_add];
%----------------------------------------------------------------------
                face_rem = zeros(6,1);
                n_face = size(mod.mod{i_mod(1)}.var.face_unq,1);
                id_all = 1:n_face;
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==i_e(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
                i_tem = sum(mod.mod{i_mod(1)}.var.face_unq==i_e(1),2) + sum(mod.mod{i_mod(1)}.var.face_unq==j(2),2) + sum(mod.mod{i_mod(1)}.var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
                id_rem = false(n_face,1);
                id_rem(face_rem) = true;
                var_new.face_unq(id_rem,:) = [];
%--------------------------------------------------------------------------
                face_add = zeros(4,3);
                i_f = 1;
                if sum(edg_add(i_f,:) == i_e(1),2) == 1
                    face_add(i_f,:) = [j(1) i_e(1) k(2)]; face_add(i_f+1,:) = [j(1) k(1) i_e(1)];
                else
                    face_add(i_f,:) = [j(1) k(1) k(2)]; face_add(i_f+1,:) = [i_e(1) k(1) k(2)];
                end
                i_f = 2;
                if sum(edg_add(i_f,:) == i_e(1),2) == 1
                    face_add(i_f*2-1,:) = [j(2) k(3) i_e(1)]; face_add(i_f*2,:) = [i_e(1) k(4) j(2)];
                else
                    face_add(i_f*2-1,:) = [j(2) k(3) k(4)]; face_add(i_f*2,:) = [k(4) k(3) i_e(1)];
                end
                
                var_new.face_unq = [var_new.face_unq;face_add];
%--------------------------------------------------------------------------                
                var_new.face_unq(var_new.face_unq==i_e(2)) = i_e(1);
                id_tem = var_new.face_unq>i_e(2);
                var_new.face_unq(id_tem) = var_new.face_unq(id_tem)-1;
                
                var_new.edge_all(var_new.edge_all==i_e(2)) = i_e(1);
                id_tem = var_new.edge_all>i_e(2);
                var_new.edge_all(id_tem) = var_new.edge_all(id_tem)-1;
                
                edg_add(edg_add>i_e(2))=edg_add(edg_add>i_e(2))-1;
            
                var_new.coord(i(2),:) = [];
                var_new.n_coord = var_new.n_coord-1;
%--------------------------------------------------------------------------                
                [mod.mod{i_mod(1)}.var] = ch.ModMembrane8AddVertex(mod.mod{i_mod(1)}.pm,mod.mod{i_mod(1)}.var,var_new.coord,var_new.edge_all,var_new.face_unq,'id_merge',i_e);
                [mod.mod{i_mod(1)}] = getUface(mod.mod{i_mod(1)});
                %[mod.mod{i_mod(1)}.var] = ch.ModMembrane8AddVertex(mod.mod{i_mod(1)}.pm,mod.mod{i_mod(1)}.var,ver_new,edge_all_new,face_new,'id_merge',i);
                edg_add = [edg_add;mod.mod{i_mod(1)}.var.edge_all(sum(mod.mod{i_mod(1)}.var.edge_all == i_e(1),2)>0,:)];
                edg_add=unique(edg_add,'row');
                [mod,relaxed] = ModMembrane8LocRelax(ch,mod,edg_add,'f_const_only',f_const_only,'local',1);
                
                mod.mod{i_mod(2)}.var.id_ModMembrane(mod.mod{i_mod(2)}.var.id_ModMembrane==i_e(2))= i_e(1);
                id_tem = mod.mod{i_mod(2)}.var.id_ModMembrane>i_e(2);
                mod.mod{i_mod(2)}.var.id_ModMembrane(id_tem)=mod.mod{i_mod(2)}.var.id_ModMembrane(id_tem)-1;
                mode_tem=mode(mod.mod{i_mod(2)}.var.id_ModMembrane);
                id_tem=(mod.mod{i_mod(2)}.var.id_ModMembrane==mode_tem)&(mod.mod{i_mod(2)}.var.occupied==false);
                id_all=(1:mod.mod{i_mod(2)}.var.n_coord)';
                id_transfer=id_all(id_tem);
                if ~isempty(id_transfer)
                id_transfer=id_transfer(1);      
                [~,ia] = setdiff(mod.mod{i_mod(1)}.var.id_on_coord,mod.mod{i_mod(2)}.var.id_ModMembrane);
                if isempty(ia)
                    error('nowhere to place adp after merging');
                else
                    r_remain=mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.id_on_coord(ia),:);
                    r_transfer=mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(2)}.var.id_ModMembrane(id_transfer),:);
                    dr2=sum((r_transfer-r_remain).^2,2);
                    [~,id_min]=min(dr2);
                    ia=ia(id_min);
                    mod.mod{i_mod(2)}.var.id_ModMembrane(id_transfer)=mod.mod{i_mod(1)}.var.id_on_coord(ia);
                end
                end
               else
                   if print_or_not==true
                   fprintf('failed: initial merging reached loop limit\n');
                   end
               end
            else
                [mod,relaxed] = ModMembrane8LocRelax(ch,mod,mod.mod{i_mod(1)}.var.edge_all(i_edg,:),'f_const_only',f_const_only,'local',1);
                if print_or_not==true
                fprintf('failed: unmergable\n');
                end
            end
%--------------------------------------------------------------------------
        end
%--------------------------------------------------------------------------
    else
        changed = false;
    end
%--------------------------------------------------------------------------    
end


