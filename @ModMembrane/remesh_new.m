function [M,remeshedAny] = remesh(m,M,varargin)
%--------------------------------------------------------------------------
        % remesh performs the remesh manipulation on @ModMembrane,
        % including splitting and merging (in manuscript), or flipping
        % input: 
        % m - @ModMembrane input
        % M - @model input
        % output:
        % remeshedAny - true: remesh done, false: remesh gave up
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
%% physical remeshing
if M.mod{M.i_mod.ModMembrane}.pm.remeshScheme==0
changed = true; 
remeshedAny=false;
while (changed == true)
    %% check for abnormal edges
    r = sqrt(sum(([M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,2),3)]...
        -[M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),1),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),2),M.mod{i_mod}.var.coord(M.mod{i_mod}.var.edge_all(:,1),3)]).^2,2));
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
    %% compute for each abnormal edge
%----------------------------------------------------------------------------------------    
    if nRemesh > 0
        iRemesh = randsample(nRemesh,1);
        i_edg = idRemesh(iRemesh);              
        if FSM(iRemesh) == 1 %means splitting
            modSave=M;
%----------------------------------------------------------------------------------------
            if split.can(i_edg)
%----------------------------------------------------------------------
               [M,relaxed] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',2);               
               if relaxed==2 %fail relax
%                 [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all(i_edg,:),'local',1);
                if print_or_not==true
                    fprintf('failed: initial extension\n');
                end
               elseif relaxed==1 %relaxed
%-------------------------------------------------------------------------- 
               j = split.j{i_edg}; k = split.k{i_edg}; id_ring_edg = split.id_ring_edg{i_edg}; edg_add_org = split.edg_add{i_edg};
               rLim=[M.TypForce.int_stored.ModMembrane.rg(1),M.TypForce.int_stored.ModMembrane.rg(end)];
               [m,remeshed,edg_add] = remeshSplitOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim); %call for splitting operation
                if remeshed==true %success
                    M.mod{i_mod}=m;
                    [M,~] = remeshLocRelax(m,M,edg_add,'local',1);
                    remeshedAny=true;
                else %fail
                   [~,id_ring_edg,~,~]=remeshRing(M.mod{i_mod},i_edg,'ring_ord', 1);
                   [M,~] = remeshLocRelax(m,M,M.mod{i_mod}.var.edge_all([i_edg;id_ring_edg],:),'local',1);
                   if print_or_not==true
                   fprintf('failed: no suitable configuration\n');
                   end
                end
%--------------------------------------------------------------------------
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
        else %merging
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
               j = merge.j{i_edg}; k = merge.k{i_edg}; id_ring_edg = merge.id_ring_edg{i_edg}; edg_add_org = merge.edg_add{i_edg};
               rLim=[M.TypForce.int_stored.ModMembrane.rg(1),M.TypForce.int_stored.ModMembrane.rg(end)];
               [m,remeshed] = remeshMergeOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim);
                
                if remeshed==true
                remeshedAny=true;
                M.mod{i_mod}=m;
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
%% CS-based remeshing
elseif M.mod{M.i_mod.ModMembrane}.pm.remeshScheme==1
    [M,remeshedAny] = remeshCSbased(m,M);
end
%==========================================================================
end