function [f,V_tot] = ModBrownianMotor_ModPolymerChain(f,M,varargin)
%--------------------------------------------------------------------------
        % ModBrownianMotor_ModPolymerChain performs the calculation of the 
        % force (@TypForce) between the two modules @ModBrownianMotor and 
        % @ModPolymerChain. mesh is used to locate near points for the
        % calculation to avoid messive search of points
        % input: 
        % f - @TypForce object
        % M - @model object including all of the modules in dynamics
        % Fname - name of forces to be included, e.g. Fname={'ModMembrane','ModMembrane_ModSubstrate'};
        % optional:
        % see variable arguments
        %   See also TimeEval
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/04
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addParameter('name', {'ModBrownianMotor';'ModPolymerChain'}, @iscell); %name of the two module evaluated
ip.parse(f,M,varargin{:});
%----------------------------------------------------------------------------------------
name=ip.Results.name;
%----------------------------------------------------------------------------------------
M.mod{M.i_mod.(name{1})}.var.idMesh = ...
ComMath.getMeshID(M.mod{M.i_mod.(name{1})}.var.coord,M.Mesh.coord, M.Mesh.range,M.Mesh.d);
M.mod{M.i_mod.(name{2})}.var.idMesh = ...
ComMath.getMeshID(M.mod{M.i_mod.(name{2})}.var.coord,M.Mesh.coord, M.Mesh.range,M.Mesh.d);
%----------------------------------------------------------------------------------------
idMeshMod1=M.mod{M.i_mod.(name{1})}.var.idMesh';

idMeshMod2 = ComMath.getMeshIDNeighbor(M.mod{M.i_mod.(name{2})}.var.idMesh',M.Mesh.range);

idMod1=[];
idMod2=[];
for i=1:27 %is nMeshFreeParticle=size(idMeshFreeParticle,1); 27 neighbors at one mesh point
IdM=(idMeshMod1-idMeshMod2(i,:)');

[idMod2Tem,idMod1Tem]=find(IdM==0);
idMod2=[idMod2;idMod2Tem];
idMod1=[idMod1;idMod1Tem];
end

f.int_comp.([name{1} '_' name{2}])=cell(2,1);
f.int_comp.([name{1} '_' name{2}]){1}=zeros(M.mod{M.i_mod.(name{1})}.var.n_coord,3);
f.int_comp.([name{1} '_' name{2}]){2}=zeros(M.mod{M.i_mod.(name{2})}.var.n_coord,3);

nPair=numel(idMod2);
Dpair=zeros(nPair,1);

for iPair=1:nPair
   Dpair(iPair)=...
   sum((M.mod{M.i_mod.(name{2})}.var.coord(idMod2(iPair),:)-M.mod{M.i_mod.(name{1})}.var.coord(idMod1(iPair),:)).^2,2);
end

[Cunq,~,iCunq]=unique(idMod1,'row');
Nunq=size(Cunq,1);
Dunq=inf(Nunq,1);
idMin=zeros(Nunq,1);
for iPair=1:nPair
    if Dunq(iCunq(iPair)) > Dpair(iPair)
        Dunq(iCunq(iPair))=Dpair(iPair);
        idMin(iCunq(iPair))=iPair;
    end
end

V_tot=0;
for iUnq=1:Nunq
    iPair=idMin(iUnq);
    [V,F] = potential.Vlinear(M.mod{M.i_mod.(name{1})}.var.coord(idMod1(iPair),:),...
                              M.mod{M.i_mod.(name{2})}.var.coord(idMod2(iPair),:),...
                             [M.mod{M.i_mod.(name{1})}.pm.A,M.mod{M.i_mod.(name{2})}.pm.r]);
    V_tot=V_tot+V;
    f.int_comp.([name{1} '_' name{2}]){1}(idMod1(iPair),:)=F;
    f.int_comp.([name{1} '_' name{2}]){2}(idMod2(iPair),:)=-F;
    if idMod2(iPair)<M.mod{M.i_mod.(name{2})}.var.n_coord
    [V,F] = potential.Vlinear(M.mod{M.i_mod.(name{1})}.var.coord(idMod1(iPair),:),...
                              M.mod{M.i_mod.(name{2})}.var.coord(idMod2(iPair)+1,:),...
                             [M.mod{M.i_mod.(name{1})}.pm.Adrag,M.mod{M.i_mod.(name{2})}.pm.r]);
    V_tot=V_tot+V;
    f.int_comp.([name{1} '_' name{2}]){1}(idMod1(iPair),:)=f.int_comp.([name{1} '_' name{2}]){1}(idMod1(iPair),:)+F;
    f.int_comp.([name{1} '_' name{2}]){2}(idMod2(iPair)+1,:)=f.int_comp.([name{1} '_' name{2}]){2}(idMod2(iPair)+1,:)-F;
    end
    if idMod2(iPair)>1
    [V,F] = potential.Vlinear(M.mod{M.i_mod.(name{1})}.var.coord(idMod1(iPair),:),...
                              M.mod{M.i_mod.(name{2})}.var.coord(idMod2(iPair)-1,:),...
                             [M.mod{M.i_mod.(name{1})}.pm.A,M.mod{M.i_mod.(name{2})}.pm.r]);
    V_tot=V_tot+V;
    f.int_comp.([name{1} '_' name{2}]){1}(idMod1(iPair),:)=f.int_comp.([name{1} '_' name{2}]){1}(idMod1(iPair),:)+F;
    f.int_comp.([name{1} '_' name{2}]){2}(idMod2(iPair)+1,:)=f.int_comp.([name{1} '_' name{2}]){2}(idMod2(iPair)-1,:)-F;
    end
end
f.int_V.([name{1} '_' name{2}])=V_tot;
f.int_tot.([name{1} '_' name{2}])=f.int_comp.([name{1} '_' name{2}]);
%--------------------------------------------------------------------------
%%
% figure();
% quiver3(mod.mod{mod.i_mod.ModMembrane}.var.coord(:,1),...
%         mod.mod{mod.i_mod.ModMembrane}.var.coord(:,2),...
%         mod.mod{mod.i_mod.ModMembrane}.var.coord(:,3),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,1),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,2),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,3),'linewidth',2);
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
%             plot(mod.mod{mod.i_mod.ModMembrane},'linestyle','-','f',fig,'FaceAlpha',1);
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,3),10,'filled');
%             hold on;
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),3),20,'filled');
%                  x_lim=[-15 15];
%            xlim(x_lim);ylim(x_lim);zlim(x_lim);
%--------------------------------------------------------------------------
end