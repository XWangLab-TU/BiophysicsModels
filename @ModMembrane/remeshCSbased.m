function [M,remeshed] = remeshCSbased(m,M,varargin)
%--------------------------------------------------------------------------
        % remeshCSbased performs the remesh manipulation on @ModMembrane,
        % including splitting and merging or flipping without free energy
        % input: 
        % m - @ModMembrane input
        % M - @model input
        % output:
        % remeshed - true: remesh done, false: remesh gave up
        % optional:
        % see variable arguments
        %   See also remeshFlip, remeshSplitMerge
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addParameter('print_or_not', true, @islogical); %whether to ouput steps
ip.parse(m,M,varargin{:});
%----------------------------------------------------------------------------------------
print_or_not=ip.Results.print_or_not;
%--------------------------------------------------------------------------
            i_mod=M.i_mod.ModMembrane;  %membrane
%----------------------------------------------------------------------------------------
    Vpm=M.mod{i_mod}.pm.Vdh;
%----------------------------------------------------------------------------------------
%%
changed = true; 
remeshed=false;
while (changed == true)
    %%
    r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
        -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
    %----------------------------------------------------------
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %----------------------------------------------------------
    idTooLong = r>Vpm.rl_max;    
    id_all = (1:M.mod{i_mod}.var.n_edg)';   
    idTooLong = id_all(idTooLong);
    

        idTooShort = r<Vpm.rl_min;
        idTooShort = id_all(idTooShort);
        [split,merge] = remeshSplitMerge(m,M,idTooLong,idTooShort);
        idRemesh = [idTooLong;idTooShort];
        FSM = [ones(length(idTooLong),1);2*ones(length(idTooShort),1)]; %Flip=0 Split=1 Merge=2
 
    
    nRemesh = size(idRemesh,1);
    if print_or_not==true
            fprintf('abnormal edge: %d, long: %d, short: %d\n',nRemesh,numel(idTooLong),numel(idTooShort));
    end
    %%
%----------------------------------------------------------------------------------------    
    if nRemesh > 0
        iRemesh = randsample(nRemesh,1);
        i_edg = idRemesh(iRemesh);
        if FSM(iRemesh) == 0
            modSave=M;
            if flip.can(i_edg)
                [M,relaxed] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',2);
                %----------------------------------------------------------
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                      -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %----------------------------------------------------------
                if relaxed==0
                   M=modSave;
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                %----------------------------------------------------------
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                      -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %----------------------------------------------------------
                   if print_or_not==true
                    fprintf('failed: initial flipping\n');
                   end
                end

                if relaxed==2
                if print_or_not==true
                    fprintf('failed: initial flipping\n');
                end
                elseif relaxed==1
                %----------------------------------------------------------------------
                remeshed=true;
                var_new=M.mod{i_mod}.var;
                j = flip.j{i_edg}; i_e = M.mod{i_mod}.var.edge_all(i_edg,:);
                %----------------------------------------------------------------------
                var_new.edge_all(i_edg,:)=[j(1) j(2)];
                if sqrt(sum((var_new.coord(j(2),:)-var_new.coord(j(1),:)).^2,2))<Vpm.r_best_min
                    disp('too short')
                end
                %----------------------------------------------------------------------
                FaceToAdd=zeros(2,3); 
                idAll=(1:size(M.mod{i_mod}.var.face_unq,1))';
                idTem=sum((M.mod{i_mod}.var.face_unq==i_e(1))+(M.mod{i_mod}.var.face_unq==j(1)),2);
                idTem=idAll(idTem==2);
                if ismember(i_e(2),M.mod{i_mod}.var.face_unq(idTem(1),:))
                    idTem(1,:)=[];
                else
                    idTem(2,:)=[];
                end
                
                nZeroMax=0;
                faceAttempt=[i_e(1) j(1) 0; 0 i_e(1) j(1); j(1) 0 i_e(1)];
%                 faceAttempt=[j(1) i_e(1) 0; 0 j(1) i_e(1); i_e(1) 0 j(1)];
                faceOrg=M.mod{i_mod}.var.face_unq(idTem,:);
                for iA=1:3
                Dface=faceOrg-faceAttempt(iA,:);
                nZero=numel(Dface(Dface==0));
                if nZeroMax<nZero
                    nZeroMax=nZero;
                end
                end
                if nZeroMax==2
                    FaceToAdd(1,:)=[i_e(1) j(2) j(1)];
                else
                    FaceToAdd(1,:)=[i_e(1) j(1) j(2)];
                end
                
                idTem=sum((M.mod{i_mod}.var.face_unq==i_e(2))+(M.mod{i_mod}.var.face_unq==j(1)),2);
                idTem=idAll(idTem==2);
                if ismember(i_e(1),M.mod{i_mod}.var.face_unq(idTem(1),:))
                    idTem(1,:)=[];
                else
                    idTem(2,:)=[];
                end
                
                nZeroMax=0;
                faceAttempt=[i_e(2) j(1) 0; 0 i_e(2) j(1); j(1) 0 i_e(2)];
                faceOrg=M.mod{i_mod}.var.face_unq(idTem,:);
                for iA=1:3
                Dface=faceOrg-faceAttempt(iA,:);
                nZero=numel(Dface(Dface==0));
                if nZeroMax<nZero
                    nZeroMax=nZero;
                end
                end
                if nZeroMax==2
                    FaceToAdd(2,:)=[i_e(2) j(2) j(1)];
                else
                    FaceToAdd(2,:)=[i_e(2) j(1) j(2)];
                end
                
                idFaceToRem=sum((M.mod{i_mod}.var.face_unq==i_e(1))+(M.mod{i_mod}.var.face_unq==i_e(2)),2);
                idFaceToRem=idAll(idFaceToRem==2);
                %----------------------------------------------------------------------
                var_new.face_unq(idFaceToRem,:)=FaceToAdd;
                %----------------------------------------------------------------------
                [M.mod{i_mod}.var,topologicalDefect] = m.remeshAddVertex(M.mod{i_mod}.pm,M.mod{i_mod}.var,var_new.coord,var_new.edge_all,var_new.face_unq);
                if topologicalDefect==true
                   M.mod{i_mod}.failInfo='topologicalDefect';
                else
                [M.mod{i_mod}] = getUface(M.mod{i_mod});
                edg_add=[j(1) j(2)];
                [M,relaxed] = remeshLocRelax(m,M,edg_add,'local',1);
                end
                %----------------------------------------------------------
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                      -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %%
%                 fig=figure('units','normalized','outerposition',[0 0 1 1]);
%                 plot(mod.mod{mod.i_mod.ModMembrane},'f',fig,'LineStyle','-');
%                 A=mod.mod{i_mod}.var.edge_all;
%                 coord=mod.mod{i_mod}.var.coord;
%                 plot3(coord(A(idMin,:),1),coord(A(idMin,:),2),coord(A(idMin,:),3),'linewidth',2);hold on;
%                 scatter3(coord(A(idMin,:),1),coord(A(idMin,:),2),coord(A(idMin,:),3),'filled')
%                 hold on;
%                 plot3(coord(edg_add,1),coord(edg_add,2),coord(edg_add,3),'linewidth',2);hold on;
                %%
                %----------------------------------------------------------
                if relaxed==0
                    M=modSave;
                    [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                     %----------------------------------------------------------
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                      -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %%
%                 [~,idMin]=min(i);
%                 edg_add=mod.mod{i_mod}.var.edge_all(i_edg,:);
%                 fig=figure('units','normalized','outerposition',[0 0 1 1]);
%                 plot(mod.mod{mod.i_mod.ModMembrane},'f',fig,'LineStyle','-');
%                 A=mod.mod{i_mod}.var.edge_all;
%                 coord=mod.mod{i_mod}.var.coord;
%                 plot3(coord(A(idMin,:),1),coord(A(idMin,:),2),coord(A(idMin,:),3),'linewidth',2);hold on;
%                 scatter3(coord(A(idMin,:),1),coord(A(idMin,:),2),coord(A(idMin,:),3),'filled')
%                 hold on;
%                 plot3(coord(edg_add,1),coord(edg_add,2),coord(edg_add,3),'linewidth',2);hold on;
                %%
                %----------------------------------------------------------
                   if print_or_not==true
                    fprintf('failed: initial flipping\n');
                   end
                else
                    if print_or_not==true
                    fprintf('successful flipping\n');
                    end
                end
                end
            else
                %----------------------------------------------------------------------
                [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                 %----------------------------------------------------------
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                      -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                i_shift=M.TypForce.int_stored.ModMembrane.rn(1)/M.mod{i_mod}.pm.dr-1;
                i = floor(r/M.mod{i_mod}.pm.dr+0.5)-i_shift;
                if min(i)<=0
                    disp(min(i));
                    disp('wall failed');
                end
                %----------------------------------------------------------
                if print_or_not==true
                fprintf('failed: unflippable\n');
                end
                %----------------------------------------------------------------------
            end           
        elseif FSM(iRemesh) == 1
            modSave=M;
%----------------------------------------------------------------------------------------
            if split.can(i_edg)
%----------------------------------------------------------------------
               [M,relaxed] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',2);               
               if relaxed==2
%                 [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                if print_or_not==true
                    fprintf('failed: initial extension\n');
                end
               elseif relaxed==1
%----------------------------------------------------------------------
                remeshed=true;
                var_new=M.mod{i_mod}.var;
                var_new.coord = [M.mod{i_mod}.var.coord;0.5*(M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(i_edg,1),:)+M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(i_edg,2),:))];
                var_new.n_coord = var_new.n_coord+1;
%----------------------------------------------------------------------
                j = split.j{i_edg}; k = split.k{i_edg}; id_ring_edg = split.id_ring_edg{i_edg}; i_n = M.mod{i_mod}.var.n_coord+1; i_e = M.mod{i_mod}.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
                dens_new=var_new.dens;
                densTem=sum(dens_new(i_e))/3;
                dens_new(i_e)=dens_new(i_e)/3*2;
                dens_new=[dens_new; densTem];
%----------------------------------------------------------------------
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
                [M.mod{i_mod}.var,topologicalDefect] = m.remeshAddVertex(M.mod{i_mod}.pm,M.mod{i_mod}.var,var_new.coord,var_new.edge_all,var_new.face_unq,...
                                                                'dens',dens_new);
                if topologicalDefect==true
                   M.mod{i_mod}.failInfo='topologicalDefect';
                else                   
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                  -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                if (max(r)<M.TypForce.int_stored.ModMembrane.rg(end)) && (min(r)>M.TypForce.int_stored.ModMembrane.rg(1))
                    successTem=true;
                    break;
                elseif (i_try==n_try)
                    successTem=false;
                end
                end
                M=modSave;
                var_new=var_new_save;
                end
                if successTem==true
                [M.mod{i_mod}] = getUface(M.mod{i_mod});
                [M,~] = remeshLocRelax(m,M,edg_add,'local',1);
                else
                   M=modSave;
                   [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all([i_edg;id_ring_edg],:),'local',1);
                   if print_or_not==true
                   fprintf('failed: initial extension too long\n');
                   end
                end
               else
                   M=modSave;
                   [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all([i_edg;id_ring_edg],:),'local',1);
                   if print_or_not==true
                   fprintf('failed: initial extension reached loop limit\n');
                   end
               end
            else
%--------------------------------------------------------------------------
                [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                if print_or_not==true
                fprintf('failed: unsplittable\n');
                end
%--------------------------------------------------------------------------
            end
%----------------------------------------------------------------------------------------
        else
            modSave=M;
            if merge.can(i_edg)
%----------------------------------------------------------------------------------------
               [M,relaxed] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',3);
               
               if relaxed==2
%                 [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                if print_or_not==true
                fprintf('failed: initial merging\n');
                end
               elseif relaxed==1
%----------------------------------------------------------------------------------------
                remeshed=true;
                var_new=M.mod{i_mod}.var;
                i = M.mod{i_mod}.var.edge_all(i_edg,:);
                var_new.coord(i(1),:) = 0.5*(M.mod{i_mod}.var.coord(i(1),:)+M.mod{i_mod}.var.coord(i(2),:));
%----------------------------------------------------------------------
                dens_new=var_new.dens;
                dens_new(i(1))=dens_new(i(1))+dens_new(i(2));
                dens_new(i(2))=[];
%----------------------------------------------------------------------
                j = merge.j{i_edg}; k = merge.k{i_edg}; id_ring_edg = merge.id_ring_edg{i_edg}; i_e = M.mod{i_mod}.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
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
                r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
                  -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
                
                if (max(r)<M.TypForce.int_stored.ModMembrane.rg(end)) && (min(r)>M.TypForce.int_stored.ModMembrane.rg(1))
                    successTem=true;
                    break;
                elseif (i_try==n_try)
                    successTem=false;
                end
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
                [M,~] = remeshLocRelax(m,M,edg_add,'local',1);
                else
                   M=modSave;
                   [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all([i_edg;id_ring_edg],:),'local',1);
                   if print_or_not==true
                   fprintf('failed: initial merging too long\n');
                   end
                end
               else
                   M=modSave;
                   [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all([i_edg;id_ring_edg],:),'local',1);
                   if print_or_not==true
                   fprintf('failed: initial merging reached loop limit\n');
                   end
               end
            else
                [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
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
