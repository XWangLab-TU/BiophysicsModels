function [M,changed] = docking(c,M,r_init,varargin)
%--------------------------------------------------------------------------
        % init performs the initialization of @ModClathrin, putting a
        % clathrin near the assigned center
        % input: 
        % c - @ModClathrin object
        % M - @model object
        % r_init - assigned center
        % optional:
        % see variable arguments
        %   See also add
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('r_init', @(x) isnumeric(x));
ip.addParameter('update', false, @islogical); %false: 1st clathrin; true: add clathrin
ip.addParameter('shiftRange', 5, @isnumeric); %defines the range of random search
ip.addParameter('nSearchXYZ', 100, @isnumeric); 
ip.addParameter('xyzLim', [-inf inf;-inf inf;-inf inf], @isnumeric);
ip.parse(c,M,r_init,varargin{:});
%----------------------------------------------------------------------------------------
i_mod=[M.i_mod.ModClathrin,M.i_mod.ModMembrane,M.i_mod.ModFreeParticle];  %1: clathrin, 2: membrane, 3: FreeParticle(AP2)
update=ip.Results.update;
shiftRange=ip.Results.shiftRange;
nSearchXYZ=ip.Results.nSearchXYZ;
xyzLim=ip.Results.xyzLim;
%----------------------------------------------------------------------------------------
%%
changed=false; %whether successfully initialized clathrin
id_feet=[49 50 51];

% distance between 1st and 2nd sphere
shift_xyz=2*(rand(nSearchXYZ,3)-0.5)*norm(M.mod{i_mod(1)}.var.coord_org(1,1:3)-M.mod{i_mod(1)}.var.coord_org(2,1:3))*shiftRange; 
% shift_xyz=zeros(1,3);

n_xyz=size(shift_xyz,1);
%--------------------------------------------------------------------------
% getting the index of vertex of @ModMembrane in the mesh defined in M
M.mod{M.i_mod.ModFreeParticle}.var.idMesh = ...
ComMath.getMeshID(M.mod{M.i_mod.ModFreeParticle}.var.coord,M.Mesh.coord, M.Mesh.range,M.Mesh.d);
%--------------------------------------------------------------------------
% random search for an acceptable location for the clathrin
initDone=false;
for i_xyz=1:n_xyz
%--------------------------------------------------------------------------
   coord_new=r_init+shift_xyz(i_xyz,:);
%--------------------------------------------------------------------------
%%
% for each random search point, rotate clathrin to add more randomness to
% find best location
for i_test=1:M.mod{i_mod(1)}.pm.n_test_all
   %-----------------------------------------------------------------------
   foot=M.mod{i_mod(1)}.get_r_from_a(M.mod{i_mod(1)}.var.coord_org(id_feet,1:3),3,M.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
   %-----------------------------------------------------------------------
   idMeshC = ComMath.getMeshID(foot,M.Mesh.coord,M.Mesh.range,M.Mesh.d);
   [Cmatch,idFP]=ismember(idMeshC,M.mod{M.i_mod.ModFreeParticle}.var.idMesh);
   Nmatch=numel(Cmatch(Cmatch==true));
   
   if Nmatch>=M.mod{i_mod(1)}.pm.n_init_adp %required number of docked feet to vertices
       %-----------------------------------------------------------------------
%        xyzLimTest=M.mod{M.i_mod.ModFreeParticle}.var.coord(idFP,:);
%        x_in=(xyzLimTest(:,1)-xyzLim(1,1)>0)&(xyzLimTest(:,1)-xyzLim(1,2)<0);
%        y_in=(xyzLimTest(:,2)-xyzLim(2,1)>0)&(xyzLimTest(:,2)-xyzLim(2,2)<0);
%        z_in=(xyzLimTest(:,3)-xyzLim(3,1)>0)&(xyzLimTest(:,3)-xyzLim(3,2)<0);
       in_or_not=true;
%        for iLimTest=1:Nmatch
%            if ((x_in(iLimTest)==true)&&(y_in(iLimTest)==true)&&(z_in(iLimTest)==true))
%            else
%                in_or_not=false;
%                break;
%            end
%        end
%--------------------------------------------------------------------------
   if in_or_not==true
       initDone=true;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
       %assign clathrin variable
       if update==false
           M.mod{i_mod(1)}.var.n_coord=1;
           M.mod{i_mod(1)}.var.connect=zeros(3,2,1);
           M.mod{i_mod(1)}.var.coord=coord_new;
           M.mod{i_mod(1)}.var.a=M.mod{i_mod(1)}.pm.a_test_all(i_test,:);
           M.mod{i_mod(1)}.var.ang_a.Phi=M.mod{i_mod(1)}.pm.ang_a_test_all(i_test,1);
           M.mod{i_mod(1)}.var.idBond=idFP';
       else
           M.mod{i_mod(1)}.var.n_coord=M.mod{i_mod(1)}.var.n_coord+1;
           M.mod{i_mod(1)}.var.connect=cat(3,M.mod{i_mod(1)}.var.connect,zeros(3,2));
           M.mod{i_mod(1)}.var.coord=[M.mod{i_mod(1)}.var.coord;coord_new];
           M.mod{i_mod(1)}.var.a=[M.mod{i_mod(1)}.var.a;M.mod{i_mod(1)}.pm.a_test_all(i_test,:)];
           M.mod{i_mod(1)}.var.ang_a.Phi=[M.mod{i_mod(1)}.var.ang_a.Phi;M.mod{i_mod(1)}.pm.ang_a_test_all(i_test,1)];
           M.mod{i_mod(1)}.var.idBond=[M.mod{i_mod(1)}.var.idBond;idFP'];
       end
       
       M.mod{i_mod(1)}.var.O=M.mod{i_mod(1)}.Omega(M.mod{i_mod(1)}.var.a,M.mod{i_mod(1)}.var.ang_a.Phi);
       M.mod{i_mod(1)}.var.E=M.mod{i_mod(1)}.getE(M.mod{i_mod(1)}.var.a,M.mod{i_mod(1)}.var.ang_a.Phi);
%--------------------------------------------------------------------------
   end
   end
   
   if initDone==true
       changed=true;
       break;
   end
end
%--------------------------------------------------------------------------
   if initDone==true
       changed=true;
       break;
   end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------    
%%
% mod.mod{mod.i_mod.ModClathrin} = ...
% getFootIDmesh(mod.mod{mod.i_mod.ModClathrin},mod.Mesh);
% idMeshClathrin=mod.mod{mod.i_mod.ModClathrin}.var.idMesh';
% IdM=(idMeshClathrin-mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh);
% 
% [idFreeParticle,idClathrin]=find(IdM==0);
% idClathrin = [floor(idClathrin/3-1e-15)+1,rem(idClathrin,3)];
% idClathrin(idClathrin(:,2)==0,2)=3;
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
% %             plot(mod.mod{mod.i_mod.ModMembrane},'linestyle','-','f',fig,'FaceAlpha',1,'col',mod.mod{mod.i_mod.ModMembrane}.var.f.kH,'col_max',0.5,'view_ang',[view_theta view_phi]);
%             col_tem=zeros(mod.mod{mod.i_mod.ModClathrin}.var.n_coord,3); col_tem(:,1)=1;col_tem(:,3)=1;
%             plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',false,'col',col_tem,...
%                 'simple',true,'iC_iLeg',idClathrin);%'iC_iLeg',[6 1;4 2] 
%             hold on;
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,3),40,'filled');
%                  x_lim=[-15 15];
%            xlim(x_lim);ylim(x_lim);zlim(x_lim);