% function [f,V_tot,identifier,idConnect] = ModClathrin_ModFreeParticle(f,mod,varargin)
function [f,V_tot,otherInfo] = ModClathrin_ModFreeParticle(f,mod,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('dr', 0.000001, @isnumeric);
ip.addParameter('idClathrinSub', [], @isnumeric);
ip.addParameter('Vonly', false, @islogical);
ip.addParameter('otherInfo', [], @isstruct);
ip.parse(f,mod,varargin{:}); 
%----------------------------------------------------------------------------------------
dr=ip.Results.dr;
Vonly=ip.Results.Vonly;
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModFreeParticle];  %1: clathrin, 2: AP2 as free particle
%----------------------------------------------------------------------------------------
idClathrinSub=ip.Results.idClathrinSub;
%----------------------------------------------------------------------------------------
%%
otherInfo=ip.Results.otherInfo;
if isempty(otherInfo)
    reload=false;
else
    reload=true;
end
%----------------------------------------------------------------------------------------
if reload==false
%----------------------------------------------------------------------------------------
mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh = ...
ComMath.getMeshID(mod.mod{mod.i_mod.ModFreeParticle}.var.coord,mod.Mesh.coord, mod.Mesh.range,mod.Mesh.d);
mod.mod{mod.i_mod.ModClathrin} = ...
getFootIDmesh(mod.mod{mod.i_mod.ModClathrin},mod.Mesh);
%----------------------------------------------------------------------------------------
if ~isempty(idClathrinSub)
    idMeshClathrinTem=mod.mod{mod.i_mod.ModClathrin}.var.idMesh';
    idMeshClathrin=-ones(size(idMeshClathrinTem));
    idMeshClathrin((idClathrinSub-1)*3+1)=idMeshClathrinTem((idClathrinSub-1)*3+1);
    idMeshClathrin((idClathrinSub-1)*3+2)=idMeshClathrinTem((idClathrinSub-1)*3+2);
    idMeshClathrin((idClathrinSub-1)*3+3)=idMeshClathrinTem((idClathrinSub-1)*3+3);
else
    idMeshClathrin=mod.mod{mod.i_mod.ModClathrin}.var.idMesh';
end

[idMeshFreeParticle] = ComMath.getMeshIDNeighbor(mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh',mod.Mesh.range);

idClathrin=[];
idFreeParticle=[];
for i=1:27 %is nMeshFreeParticle=size(idMeshFreeParticle,1); 27 neighbors at one mesh point
IdM=(idMeshClathrin-idMeshFreeParticle(i,:)');

[idFreeParticleTem,idClathrinTem]=find(IdM==0);
idClathrinTem = [floor(idClathrinTem/3-0.1)+1,rem(idClathrinTem,3)];
idClathrinTem(idClathrinTem(:,2)==0,2)=3;
idClathrin=[idClathrin;idClathrinTem];
idFreeParticle=[idFreeParticle;idFreeParticleTem];
end
identifier=zeros(mod.n_mod,1);
identifier(mod.i_mod.ModClathrin)=1;
identifier(mod.i_mod.ModFreeParticle)=2;
otherInfo=struct('idClathrin',idClathrin,'idFreeParticle',idFreeParticle,'identifier',identifier);
else
idClathrin=otherInfo.idClathrin;
idFreeParticle=otherInfo.idFreeParticle;
end

f.int_comp.ModClathrin_ModFreeParticle=cell(2,1);
f.int_comp.ModClathrin_ModFreeParticle{1}=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
f.int_comp.ModClathrin_ModFreeParticle{2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);

nPair=numel(idFreeParticle);
coordFoot=getFoot(mod.mod{mod.i_mod.ModClathrin});
Vpair=zeros(nPair,1);
Fpair=zeros(nPair,3);

for iPair=1:nPair
   Vpair(iPair)=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((coordFoot(idClathrin(iPair,2),:,idClathrin(iPair,1))-mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)).^2,2);
   Fpair(iPair,:) = -f.pm.k_ModClathrin_ModFreeParticle*(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)-...
                                               coordFoot(idClathrin(iPair,2),:,idClathrin(iPair,1)));
%    f.int_comp.ModClathrin_ModFreeParticle{2}(idFreeParticle(iPair),:)=f.int_comp.ModClathrin_ModFreeParticle{2}(idFreeParticle(iPair),:)-f_tem;
end

[Cunq,~,iCunq]=unique(idClathrin,'row');

Nunq=size(Cunq,1);

Vunq=inf(Nunq,1);

idMin=zeros(Nunq,1);

for iPair=1:nPair
    if Vunq(iCunq(iPair)) > Vpair(iPair)
        Vunq(iCunq(iPair))=Vpair(iPair);
        idMin(iCunq(iPair))=iPair;
    end
end
for iUnq=1:Nunq
    iPair=idMin(iUnq);
    f.int_comp.ModClathrin_ModFreeParticle{2}(idFreeParticle(iPair),:)=Fpair(iPair,:);
end

% idConnect=[idClathrin(idMin,:),idFreeParticle(idMin)];
V_tot=sum(Vunq);
f.int_V.ModClathrin_ModFreeParticle=V_tot;
%%
if Vonly==false
id_feet=[49 50 51];
for iUnq=1:Nunq
    iPair=idMin(iUnq);
    for m=1:3
%--------------------------------------------------------------------------
       coordFootTem=coordFoot; 
       Vtem=V_tot-Vpair(iPair);
       coordFootTem(idClathrin(iPair,2),m,idClathrin(iPair,1))=coordFootTem(idClathrin(iPair,2),m,idClathrin(iPair,1))+dr;
%--------------------------------------------------------------------------   
       VpairTem=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((coordFootTem(idClathrin(iPair,2),:,idClathrin(iPair,1))-mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)).^2,2);
       Vtem=Vtem+VpairTem;
       f.int_comp.ModClathrin_ModFreeParticle{1}(idClathrin(iPair,1),m)=...
              f.int_comp.ModClathrin_ModFreeParticle{1}(idClathrin(iPair,1),m)-(Vtem-V_tot)/dr;
%--------------------------------------------------------------------------
       c=mod.mod{i_mod(1)};
       Vtem=V_tot-Vpair(iPair);
%--------------------------------------------------------------------------       
       c.var.a(idClathrin(iPair,1),m)=c.var.a(idClathrin(iPair,1),m)+dr;
       c.var.ang_a.Phi(idClathrin(iPair,1))=norm(c.var.a(idClathrin(iPair,1),:));
       [c.var] = setVar(c,true);
       foot_tem=c.get_r_from_a(c.var.coord_org(id_feet(idClathrin(iPair,2)),1:3),1,c.var.O(:,:,idClathrin(iPair,1)))+c.var.coord(idClathrin(iPair,1),:);    
       VpairTem=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((foot_tem-mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)).^2,2);
       Vtem=Vtem+VpairTem;
       f.int_comp.ModClathrin_ModFreeParticle{1}(idClathrin(iPair,1),m+3)=...
              f.int_comp.ModClathrin_ModFreeParticle{1}(idClathrin(iPair,1),m+3)-(Vtem-V_tot)/dr;
%--------------------------------------------------------------------------
    end
end
f.int_tot.ModClathrin_ModFreeParticle=f.int_comp.ModClathrin_ModFreeParticle;
end
%--------------------------------------------------------------------------
%%
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
% %             plot(mod.mod{mod.i_mod.ModMembrane},'linestyle','-','f',fig,'FaceAlpha',1,'col',mod.mod{mod.i_mod.ModMembrane}.var.f.kH,'col_max',0.5,'view_ang',[view_theta view_phi]);
%             col_tem=zeros(mod.mod{mod.i_mod.ModClathrin}.var.n_coord,3); col_tem(:,1)=1;col_tem(:,3)=1;
%             plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',false,'col',col_tem,...
%                 'simple',true,'iC_iLeg',idClathrin);%'iC_iLeg',[6 1;4 2] 
%             hold on;
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,3),40,'filled');
%             hold on;
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),3),60,'filled');
%                  x_lim=[-15 15];
%            xlim(x_lim);ylim(x_lim);zlim(x_lim);
%--------------------------------------------------------------------------
end


