function [M,remeshed] = remeshLattice(m,lc,M,varargin)
%--------------------------------------------------------------------------
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2025/08/27
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('lc', @(x) isa(x,'ComLattice'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addParameter('print_or_not', true, @islogical); %whether to ouput steps
ip.addParameter('nTryMax', 5, @isnumeric);
ip.parse(m,lc,M,varargin{:});
%----------------------------------------------------------------------------------------
nTryMax=ip.Results.nTryMax;
%--------------------------------------------------------------------------
i_mod=M.i_mod.ModMembrane;  %membrane
%----------------------------------------------------------------------------------------        
%========================================================================================
nTry=0;
changed=true;
while (changed == true)
%----------------------------------------------------------------------------------------
        [idRemesh,SplitOrMerge,split,merge]=getIDremesh(M,lc,i_mod);
        nRemesh = size(idRemesh,1);
%         disp(nRemesh);
%----------------------------------------------------------------------------------------        
        if nRemesh > 0
            iRemesh = randsample(nRemesh,1);
            i_edg = idRemesh(iRemesh);
            %--------------------------------------------------------------
            if SplitOrMerge(iRemesh) == 1 %split
                SplitOrMergePair=[1 2];
            else
                SplitOrMergePair=[2 1];
            end
            for iSM=1:2
                if iSM==2
                    [i_edg,split,merge]=getIDremeshPaired(M,lc,i_mod,SplitOrMergePair(iSM));
                end
                switch SplitOrMergePair(iSM)
%----------------------------------------------------------------------------------------                    
                    case 1 %split
                    modSave=M;
                    %------------------------------------------------------
                    var_new=M.mod{i_mod}.var;
                    var_new.coord = [M.mod{i_mod}.var.coord;0.5*(M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(i_edg,1),:)+M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(i_edg,2),:))];
                    var_new.n_coord = var_new.n_coord+1;
                    j = split.j{i_edg}; k = split.k{i_edg}; id_ring_edg = split.id_ring_edg{i_edg}; i_n = M.mod{i_mod}.var.n_coord+1; i_e = M.mod{i_mod}.var.edge_all(i_edg,:);
                    dens_new=var_new.dens;
                    densTem=sum(dens_new(i_e))/3;
                    dens_new(i_e)=dens_new(i_e)/3*2;
                    dens_new=[dens_new; densTem];
                    edg_rem = [i_edg;id_ring_edg(sum(M.mod{i_mod}.var.edge_all(id_ring_edg,1)==j(1),2)+sum(M.mod{i_mod}.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
                        id_ring_edg(sum(M.mod{i_mod}.var.edge_all(id_ring_edg,1)==j(2),2)+sum(M.mod{i_mod}.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
                    id_rem = false(M.mod{i_mod}.var.n_edg,1);
                    id_rem(edg_rem) = true;
                    var_new.edge_all(id_rem,:) = [];
                    edg_tem = split.edg_add{i_edg};
                    r_tem = [sqrt(sum((var_new.coord(edg_tem(:,1),:)-var_new.coord(edg_tem(:,2),:)).^2,2)),...
                        sqrt(sum((var_new.coord(edg_tem(:,3),:)-var_new.coord(edg_tem(:,4),:)).^2,2)),...
                        sqrt(sum((var_new.coord(edg_tem(:,5),:)-var_new.coord(edg_tem(:,6),:)).^2,2)),...
                        sqrt(sum((var_new.coord(edg_tem(:,7),:)-var_new.coord(edg_tem(:,8),:)).^2,2))];
                    [~,id_sort] = sort(std(r_tem-M.mod{i_mod}.pm.l0,[],2));
                    n_try=numel(id_sort);
                    successTem=false;
                    edg_save=edg_tem;
                    var_new_save=var_new;
                    for i_try=1:n_try
                        id_tem=id_sort(i_try);
                        edg_tem = edg_save(id_tem,:);
                        edg_tem = [edg_tem(1:2);edg_tem(3:4);edg_tem(5:6);edg_tem(7:8)];
                        edg_add = [edg_tem;[M.mod{i_mod}.var.edge_all(i_edg,1) i_n];[i_n M.mod{i_mod}.var.edge_all(i_edg,2)];[i_n j(1)];[i_n j(2)]];
                        var_new.edge_all = [var_new.edge_all;edg_add];
                        face_rem = zeros(6,1);
                        n_face = size(M.mod{i_mod}.var.face_unq,1);
                        id_all = 1:n_face;
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2) + sum(M.mod{i_mod}.var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==k(2),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2) + sum(M.mod{i_mod}.var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2) + sum(M.mod{i_mod}.var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
                        id_rem = false(n_face,1);
                        id_rem(face_rem) = true;
                        var_new.face_unq(id_rem,:) = [];
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
                        [M.mod{i_mod}.var,topologicalDefect] = m.remeshAddVertex(M.mod{i_mod}.pm,M.mod{i_mod}.var,var_new.coord,var_new.edge_all,var_new.face_unq,...
                            'dens',dens_new);
                        if topologicalDefect==true
                            M.mod{i_mod}.failInfo='topologicalDefect';
                        else
                            successTem=true;
                            break;
                            % r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                            %   -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                            % if (max(r)<M.TypForce.int_stored.ModMembrane.rg(end)) && (min(r)>M.TypForce.int_stored.ModMembrane.rg(1))
                            %     successTem=true;
                            %     break;
                            % elseif (i_try==n_try)
                            %     successTem=false;
                            % end
                        end
                        M=modSave;
                        var_new=var_new_save;
                    end
                    %----------------------------------------------------------
                    if successTem==true
                        [M.mod{i_mod}] = getUface(M.mod{i_mod});
                        lc = coordToMesh(lc,M);
                        [M] = remeshLocRelaxLattice(M.mod{i_mod},M,lc,edg_add);
                        lc = coordToMesh(lc,M);
                    else 
                        warning('remesh topologicalDefect');
                    end
%----------------------------------------------------------------------------------------                    
                    case 2 %merge
                    modSave=M;
                    %------------------------------------------------------
                    var_new=M.mod{i_mod}.var;
                    i = M.mod{i_mod}.var.edge_all(i_edg,:);
                    var_new.coord(i(1),:) = 0.5*(M.mod{i_mod}.var.coord(i(1),:)+M.mod{i_mod}.var.coord(i(2),:));
                    dens_new=var_new.dens;
                    dens_new(i(1))=dens_new(i(1))+dens_new(i(2));
                    dens_new(i(2))=[];
                    j = merge.j{i_edg}; k = merge.k{i_edg}; id_ring_edg = merge.id_ring_edg{i_edg}; i_e = M.mod{i_mod}.var.edge_all(i_edg,:);
                    edg_rem = [i_edg;...
                        id_ring_edg(sum(M.mod{i_mod}.var.edge_all(id_ring_edg,1)==j(1),2)+sum(M.mod{i_mod}.var.edge_all(id_ring_edg,2)==j(1),2) > 0);...
                        id_ring_edg(sum(M.mod{i_mod}.var.edge_all(id_ring_edg,1)==j(2),2)+sum(M.mod{i_mod}.var.edge_all(id_ring_edg,2)==j(2),2) > 0)];
                    %                 edg_replace = id_ring_edg(sum(M.mod{i_mod}.var.edge_all(id_ring_edg,1)==i_e(2),2)+sum(M.mod{i_mod}.var.edge_all(id_ring_edg,2)==i_e(2),2) > 0);
                    %                 [~,id_duplicated]=intersect(edg_replace,edg_rem);
                    %                 edg_replace(id_duplicated)=[];
                    %                 edg_tem=var_new.edge_all(edg_replace,:);
                    %                 edg_tem(edg_tem==i_e(2))=i_e(1);
                    %                 var_new.edge_all(edg_replace,:)=edg_tem;
                    id_rem = false(M.mod{i_mod}.var.n_edg,1);
                    id_rem(edg_rem) = true;
                    var_new.edge_all(id_rem,:) = [];
                    edg_tem = merge.edg_add{i_edg};
                    r_tem = [sqrt(sum((var_new.coord(edg_tem(:,1),:)-var_new.coord(edg_tem(:,2),:)).^2,2)),...
                        sqrt(sum((var_new.coord(edg_tem(:,3),:)-var_new.coord(edg_tem(:,4),:)).^2,2))];
                    [~,id_sort] = sort(std(r_tem-M.mod{i_mod}.pm.l0,[],2));
                    n_try=numel(id_sort);
                    successTem=false;
                    edg_save=edg_tem;
                    var_new_save=var_new;
                    for i_try=1:n_try
                        id_tem=id_sort(i_try);
                        edg_tem = edg_save(id_tem,:);
                        edg_tem = [edg_tem(1:2);edg_tem(3:4)];
                        edg_add = edg_tem;
                        var_new.edge_all = [var_new.edge_all;edg_add];
                        %----------------------------------------------------------------------
                        face_rem = zeros(6,1);
                        n_face = size(M.mod{i_mod}.var.face_unq,1);
                        id_all = 1:n_face;
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2) + sum(M.mod{i_mod}.var.face_unq==k(1),2); face_rem(1) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2); face_rem(2) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==k(2),2) + sum(M.mod{i_mod}.var.face_unq==j(1),2); face_rem(3) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2) + sum(M.mod{i_mod}.var.face_unq==k(3),2); face_rem(4) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==i_e(2),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2); face_rem(5) = id_all(i_tem==3);
                        i_tem = sum(M.mod{i_mod}.var.face_unq==i_e(1),2) + sum(M.mod{i_mod}.var.face_unq==j(2),2) + sum(M.mod{i_mod}.var.face_unq==k(4),2); face_rem(6) = id_all(i_tem==3);
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

                        [M.mod{i_mod}.var,topologicalDefect] = m.remeshAddVertex(M.mod{i_mod}.pm,M.mod{i_mod}.var,var_new.coord,var_new.edge_all,var_new.face_unq,'dens',dens_new);
                        if topologicalDefect==true
                            M.mod{i_mod}.failInfo='topologicalDefect';
                        else
                            successTem=true;
                            break;
                            % r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                            %     -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                            % 
                            % if (max(r)<M.TypForce.int_stored.ModMembrane.rg(end)) && (min(r)>M.TypForce.int_stored.ModMembrane.rg(1))
                            %     successTem=true;
                            %     break;
                            % elseif (i_try==n_try)
                            %     successTem=false;
                            % end
                        end
                        M=modSave;
                        var_new=var_new_save;
                    end
                    if successTem==true
                        [M.mod{i_mod}] = getUface(M.mod{i_mod});
                        edg_add = M.mod{i_mod}.var.edge_all(end-1:end,:);
                        i_e_new=edg_add(1,1);
                        if (edg_add(2,1)~=i_e_new) && (edg_add(2,2)~=i_e_new)
                            i_e_new=edg_add(1,2);
                        end
                        edg_add=[edg_add;M.mod{i_mod}.var.edge_all(sum(M.mod{i_mod}.var.edge_all == i_e_new,2)>0,:)];
                        edg_add=unique(edg_add,'row');
                        lc = coordToMesh(lc,M);
                        [M] = remeshLocRelaxLattice(M.mod{i_mod},M,lc,edg_add);
                        lc = coordToMesh(lc,M);
                    end
                    %------------------------------------------------------
                end
                %------------------------------------------------------
            end
%----------------------------------------------------------------------------------   
        else
            changed=false;
        end % nRemesh > 0
%----------------------------------------------------------------------------------cutoff
    nTry=nTry+1;
    disp(['abnormal edge #: ' num2str(nRemesh)]);
    if nTry>nTryMax
    warning('remesh too many trials!!');
    changed=false;
    end
%----------------------------------------------------------------------------------------    
end
%==========================================================================
end
%==========================================================================
function [idRemesh,SplitOrMerge,split,merge]=getIDremesh(M,lc,i_mod)
        meshIDedg=lc.component{i_mod}.meshID(M.mod{i_mod}.var.edge_all);
        coordTem1=meshToCoord(lc,meshIDedg(:,1));
        coordTem2=meshToCoord(lc,meshIDedg(:,2));
        dist=sqrt(sum((coordTem1-coordTem2).^2,2));
        id_all = (1:M.mod{i_mod}.var.n_edg)';   
        Rremesh1=0.5*(M.mod{i_mod}.pm.Vedg.rb_1+M.mod{i_mod}.pm.Vedg.r_1);
        Rremesh2=0.5*(M.mod{i_mod}.pm.Vedg.rb_2+M.mod{i_mod}.pm.Vedg.r_2);
        idTooShort=dist<Rremesh1;
        idTooShort=id_all(idTooShort);
        idTooLong=dist>Rremesh2;
        idTooLong=id_all(idTooLong);
        [split,merge] = remeshSplitMerge(M.mod{i_mod},M,idTooLong,idTooShort);
        %M.mod{i_mod}.var.CanSplitMerge=[split.can,merge.can];
        idTooLong=idTooLong(split.can(idTooLong));
        idTooShort=idTooShort(merge.can(idTooShort));
        idRemesh = [idTooLong;idTooShort];
        nTooshort=length(idTooShort);
        nTooLong=length(idTooLong);
        SplitOrMerge = [ones(nTooLong,1);2*ones(nTooshort,1)]; % Split=1 Merge=2
end
%==========================================================================
function [i_edg,split,merge]=getIDremeshPaired(M,lc,i_mod,minORmax)
        meshIDedg=lc.component{i_mod}.meshID(M.mod{i_mod}.var.edge_all);
        coordTem1=meshToCoord(lc,meshIDedg(:,1));
        coordTem2=meshToCoord(lc,meshIDedg(:,2));
        dist=sqrt(sum((coordTem1-coordTem2).^2,2));
        id_all = (1:M.mod{i_mod}.var.n_edg)';
        [split,merge] = remeshSplitMerge(M.mod{i_mod},M,id_all,id_all);
        if minORmax==2 %merge: find shortest
            [~,idSrt]=sort(dist,'ascend');
            idSrt=idSrt(merge.can(idSrt));
        elseif minORmax==1 %split: find shortest
            [~,idSrt]=sort(dist,'descend');
            idSrt=idSrt(split.can(idSrt));
        end
        i_edg=idSrt(1);
end
%==========================================================================