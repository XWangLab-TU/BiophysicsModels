function [M,changed] = addNear(c,M,r_init,varargin)
%--------------------------------------------------------------------------
        % init performs the addition of one more clathrin in @ModClathrin, 
        % putting one more clathrin near a leg of existing clathrin subject
        % also to an assigned center 
        % input: 
        % c - @ModClathrin object
        % M - @model object
        % r_init - the assigned center
        % optional:
        % see variable arguments
        %   See also init
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('r_init', @(x) isnumeric(x));
ip.addParameter('dt', 0.01, @isnumeric);
ip.addParameter('shiftRange', 5, @isnumeric); %defines the range of random search
ip.addParameter('nMaxRing', 6, @isnumeric); %defines the maximal number of clathrin that can be connected as a ring
ip.addParameter('sortFromCtr', true, @islogical); %whether start searching from the assigned center
ip.addParameter('rMaxFromCtr', 5, @isnumeric); %maximal distance allowed from the assigned center
ip.addParameter('nSearchXYZ', 100, @isnumeric);
ip.addParameter('xyzLim', [-inf inf;-inf inf;-inf inf], @isnumeric);
ip.parse(c,M,r_init,varargin{:});
%----------------------------------------------------------------------------------------
i_mod=[M.i_mod.ModClathrin,M.i_mod.ModMembrane,M.i_mod.ModFreeParticle];  %1: clathrin, 2: membrane 3: FreeParticle(AP2)
dt=ip.Results.dt;
shiftRange=ip.Results.shiftRange;
nMaxRing=ip.Results.nMaxRing;
sortFromCtr=ip.Results.sortFromCtr;
rMaxFromCtr=ip.Results.rMaxFromCtr;
nSearchXYZ=ip.Results.nSearchXYZ;
xyzLim=ip.Results.xyzLim;
%----------------------------------------------------------------------------------------
%%
%%
[k_on_list,n_list]=getKonList(M,i_mod,dt,sortFromCtr,r_init,rMaxFromCtr);
%%
changed=false; %whether successfully added clathrin
id_waist=[8 16 24];
id_feet=[49 50 51];

% distance between 1st and 2nd sphere
shift_xyz=2*(rand(nSearchXYZ,3)-0.5)*norm(M.mod{i_mod(1)}.var.coord_org(1,1:3)-M.mod{i_mod(1)}.var.coord_org(2,1:3))*shiftRange;

n_xyz=size(shift_xyz,1);
%--------------------------------------------------------------------------
% getting the index of vertex of @ModMembrane in the mesh defined in M
M.mod{M.i_mod.ModMembrane}.var.idMesh = ...
ComMath.getMeshID(M.mod{M.i_mod.ModMembrane}.var.coord,M.Mesh.coord, M.Mesh.range,M.Mesh.d);
%--------------------------------------------------------------------------
addDone=false;
% searching for best location of additional clathrin
for i_list=1:n_list
    
    i_c=k_on_list(i_list,1);
    i_leg_host=k_on_list(i_list,2);
%--------------------------------------------------------------------------
    junction=[i_c i_leg_host];
    % getting topology, cannot have too many clathrin connecting as a
    % circle
    [~,n_ring] = getTopology(M.mod{i_mod(1)},junction);
    if (n_ring(1) < nMaxRing) && (n_ring(2) < nMaxRing)
        waist_host=M.mod{i_mod(1)}.get_r_from_a(M.mod{i_mod(1)}.var.coord_org(id_waist(i_leg_host),1:3),1,M.mod{i_mod(1)}.var.O(:,:,i_c))+M.mod{i_mod(1)}.var.coord(i_c,:);
%--------------------------------------------------------------------------       
      for i_xyz=1:n_xyz
          coord_new=waist_host+shift_xyz(i_xyz,:);
          coord_host_hub=M.mod{i_mod(1)}.var.coord_org(end,1:3)+M.mod{i_mod(1)}.var.coord(i_c,:);
          for i_test=1:M.mod{i_mod(1)}.pm.n_test_all
%--------------------------------------------------------------------------
              foot=M.mod{i_mod(1)}.get_r_from_a(M.mod{i_mod(1)}.var.coord_org(id_feet,1:3),3,M.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
              waist=M.mod{i_mod(1)}.get_r_from_a(M.mod{i_mod(1)}.var.coord_org(id_waist,1:3),3,M.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;

              idMeshC = ComMath.getMeshID(foot,M.Mesh.coord,M.Mesh.range,M.Mesh.d);
              [Cmatch,idMem]=ismember(idMeshC,M.mod{M.i_mod.ModMembrane}.var.idMesh);
              Nmatch=numel(Cmatch(Cmatch==true));
              
              if (Nmatch>=M.mod{i_mod(1)}.pm.n_min_adp) %required number of docked feet to vertices 
                  %--------------------------------------------------------------------------
                  xyzLimTest=M.mod{M.i_mod.ModMembrane}.var.coord(idMem,:);
                  x_in=(xyzLimTest(:,1)-xyzLim(1,1)>0)&(xyzLimTest(:,1)-xyzLim(1,2)<0);
                  y_in=(xyzLimTest(:,2)-xyzLim(2,1)>0)&(xyzLimTest(:,2)-xyzLim(2,2)<0);
                  z_in=(xyzLimTest(:,3)-xyzLim(3,1)>0)&(xyzLimTest(:,3)-xyzLim(3,2)<0);
                  in_or_not=true;
                  for iLimTest=1:Nmatch
                      if ((x_in(iLimTest)==true)&&(y_in(iLimTest)==true)&&(z_in(iLimTest)==true))
                      else
                          in_or_not=false;
                          break;
                      end
                  end
              %--------------------------------------------------------------------------
              if in_or_not==true
                  addDone=true;
              %------------------------------------------------------------
                  %assign new AP2 particles to @ModFreeParticle using matched index,
                  %match ModFreeParticle to ModMembrane's vertex positions
                  M.mod{M.i_mod.ModFreeParticle}.var.n_coord=M.mod{M.i_mod.ModFreeParticle}.var.n_coord+Nmatch;
                  M.mod{M.i_mod.ModFreeParticle}.var.coord=[M.mod{M.i_mod.ModFreeParticle}.var.coord;...
                                                            M.mod{M.i_mod.ModMembrane}.var.coord(idMem,:)];
                  %--------------------------------------------------------
                  M.mod{i_mod(1)}.var.n_coord=M.mod{i_mod(1)}.var.n_coord+1;
                  M.mod{i_mod(1)}.var.connect=cat(3,M.mod{i_mod(1)}.var.connect,zeros(3,2));
                  M.mod{i_mod(1)}.var.coord=[M.mod{i_mod(1)}.var.coord;coord_new];
                  M.mod{i_mod(1)}.var.a=[M.mod{i_mod(1)}.var.a;M.mod{i_mod(1)}.pm.a_test_all(i_test,:)];
                  M.mod{i_mod(1)}.var.ang_a.Phi=[M.mod{i_mod(1)}.var.ang_a.Phi;M.mod{i_mod(1)}.pm.ang_a_test_all(i_test,1)];
                  M.mod{i_mod(1)}.var.O=M.mod{i_mod(1)}.Omega(M.mod{i_mod(1)}.var.a,M.mod{i_mod(1)}.var.ang_a.Phi);
                  M.mod{i_mod(1)}.var.E=M.mod{i_mod(1)}.getE(M.mod{i_mod(1)}.var.a,M.mod{i_mod(1)}.var.ang_a.Phi);
%--------------------------------------------------------------------------
                  M.mod{i_mod(1)}.var.connect(i_leg_host,1,i_c)=M.mod{i_mod(1)}.var.n_coord;       
                  [~,i_leg_all]=min(vecnorm(waist-coord_host_hub,2,2));
                  M.mod{i_mod(1)}.var.connect(i_leg_host,2,i_c)=i_leg_all;
                  M.mod{i_mod(1)}.var.connect(i_leg_all,1,M.mod{i_mod(1)}.var.n_coord)=i_c;
                  M.mod{i_mod(1)}.var.connect(i_leg_all,2,M.mod{i_mod(1)}.var.n_coord)=i_leg_host;
              end
              end

              if addDone==true
                  changed=true;
%                   plot(M.mod{i_mod(1)});
%                   M.mod{M.i_mod.ModFreeParticle}.var.coord
%                   pause;
                  break;
              end
          end
          if addDone==true
              changed=true;
              break;
          end
      end
    end
    if addDone==true
        changed=true;
        break;
    end
end
%--------------------------------------------------------------------------    
%%
end
function [k_on_list,n_list]=getKonList(M,i_mod,dt,sortFromCtr,r_init,rMaxFromCtr)
    k_on_list=[];
    id_waist=[8 16 24];
    r=[];
    for i_c=1:M.mod{i_mod(1)}.var.n_coord
        for i_leg=1:3
            waist=M.mod{i_mod(1)}.get_r_from_a(M.mod{i_mod(1)}.var.coord_org(id_waist(i_leg),1:3),1,M.mod{i_mod(1)}.var.O(:,:,i_c))...
                       +M.mod{i_mod(1)}.var.coord(i_c,:);
            rTem=norm(waist-r_init);
            if (rTem < rMaxFromCtr)
            r=[r;rTem]; 
            k_on_list=[k_on_list;[i_c i_leg]];
            end
        end
    end
    n_list=size(k_on_list,1);
    if sortFromCtr==true
        [~,idSortR]=sort(r);
        k_on_list=k_on_list(idSortR,:);
    else
        idRand=randsample(n_list,n_list);
        k_on_list=k_on_list(idRand,:);
    end
    idCanOn=false(n_list,1);
    for i=1:n_list
        i_c=k_on_list(i,1);
        i_leg=k_on_list(i,2);
        if (rand(1,1) < dt*M.mod{i_mod(1)}.pm.k_on) && (M.mod{i_mod(1)}.var.connect(i_leg,1,i_c) == 0) 
            idCanOn(i)=true;
        end
    end
    k_on_list=k_on_list(idCanOn,:);
    n_list=size(k_on_list,1);
end