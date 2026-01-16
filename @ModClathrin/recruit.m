function [M,changed] = recruit(c,M,varargin)
%--------------------------------------------------------------------------
        % recruit performs the addition of one more clathrin in @ModClathrin, 
        % putting one more clathrin near a leg of existing clathrin subject
        % input: 
        % c - @ModClathrin object
        % M - @model object
        % Author: Xinxin Wang
        % email: wangxinxin8627@gmail.com
        % date: 2026/01/11
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addParameter('dt', 0.01, @isnumeric);
ip.addParameter('nMaxRing', 6, @isnumeric); %defines the maximal number of clathrin that can be connected as a ring
ip.addParameter('sortFromCtr', true, @islogical); %whether start searching from the assigned center
ip.addParameter('rMaxFromCtr', 5, @isnumeric); %maximal distance allowed from the assigned center
ip.addParameter('rMinFromCtr', 0, @isnumeric); %minimal distance allowed from the assigned center
ip.addParameter('nSearchXYZ', 100, @isnumeric);
ip.addParameter('xyzLim', [-inf inf;-inf inf;-inf inf], @isnumeric);
ip.parse(c,M,varargin{:});
%----------------------------------------------------------------------------------------
[waist] = getWaist(M.mod{M.i_mod.ModClathrin});
r_init=mean(M.mod{M.i_mod.ModClathrin}.var.coord,1);
[k_on_list,n_list]=getKonList(M,ip.Results.sortFromCtr,r_init);
%----------------------------------------------------------------------------------------
addDone=false;
for i_list=1:n_list
    i_c_host=k_on_list(i_list,1);
    i_leg_host=k_on_list(i_list,2);
    junction=[i_c_host i_leg_host];
    % getting topology, cannot have too many clathrin connecting in a circle
    [~,n_ring] = getTopology(M.mod{M.i_mod.ModClathrin},junction);
    if (n_ring(1) < ip.Results.nMaxRing) && (n_ring(2) < ip.Results.nMaxRing)
        [M,changed] = docking(M.mod{M.i_mod.ModClathrin},M,waist(i_leg_host,:,i_c_host),'shiftRange',1,'update', true);
        i_c_guest=M.mod{M.i_mod.ModClathrin}.var.n_coord;
        M.mod{M.i_mod.ModClathrin}.var.connect(i_leg_host,1,i_c_host)=i_c_guest;
        [waist] = getWaist(M.mod{M.i_mod.ModClathrin});
        r=vecnorm(M.mod{M.i_mod.ModClathrin}.var.coord(i_c_host,:)-waist(:,:,i_c_guest),2,2);
        [~,i_leg_guest]=min(r);
        M.mod{M.i_mod.ModClathrin}.var.connect(i_leg_host,2,i_c_host)=i_leg_guest;
        M.mod{M.i_mod.ModClathrin}.var.connect(i_leg_guest,1,i_c_guest)=i_c_host;
        M.mod{M.i_mod.ModClathrin}.var.connect(i_leg_guest,2,i_c_guest)=i_leg_host;
        addDone=true;
    end
    if addDone==true
        changed=true;
        break;
    end
end
%----------------------------------------------------------------------------------------
end
%----------------------------------------------------------------------------------------
function [k_on_list,n_list]=getKonList(M,sortFromCtr,r_init)
    k_on_list=[];
    id_waist=[8 16 24];
    r=[];
    for i_c=1:M.mod{M.i_mod.ModClathrin}.var.n_coord
        for i_leg=1:3
            waist=M.mod{M.i_mod.ModClathrin}.get_r_from_a(M.mod{M.i_mod.ModClathrin}.var.coord_org(id_waist(i_leg),1:3),1,M.mod{M.i_mod.ModClathrin}.var.O(:,:,i_c))...
                       +M.mod{M.i_mod.ModClathrin}.var.coord(i_c,:);
            rTem=norm(waist-r_init);
            r=[r;rTem]; 
            k_on_list=[k_on_list;[i_c i_leg]];
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
        % if (rand(1,1) < dt*M.mod{M.i_mod.ModClathrin}.pm.k_on) && (M.mod{M.i_mod.ModClathrin}.var.connect(i_leg,1,i_c) == 0) 
        if (M.mod{M.i_mod.ModClathrin}.var.connect(i_leg,1,i_c) == 0) 
            idCanOn(i)=true;
        end
    end
    k_on_list=k_on_list(idCanOn,:);
    n_list=size(k_on_list,1);
end