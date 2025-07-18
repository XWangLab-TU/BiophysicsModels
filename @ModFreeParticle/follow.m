function [p,idConnect] = follow(p,mod,name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('p', @(x) isa(x,'ModFreeParticle'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('name', @(x) ischar(x));
ip.parse(p,mod,name,varargin{:}); 
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModFreeParticle,mod.i_mod.(name)];  %1: free particle, 2: mod to follow
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%%
mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh = ...
ComMath.getMeshID(mod.mod{mod.i_mod.ModFreeParticle}.var.coord,mod.Mesh.coord, mod.Mesh.range,mod.Mesh.d);
mod.mod{i_mod(2)}.var.idMesh = ...
ComMath.getMeshID(mod.mod{i_mod(2)}.var.coord,mod.Mesh.coord, mod.Mesh.range,mod.Mesh.d);
%----------------------------------------------------------------------------------------
idMeshFreeParticle=mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh';
idMeshToFollow=ComMath.getMeshIDNeighbor(mod.mod{i_mod(2)}.var.idMesh',mod.Mesh.range);

idFreeParticle=[];
idToFollow=[];
for i=1:27
IdM=(idMeshFreeParticle-idMeshToFollow(i,:)');

[idToFollowTem,idFreeParticleTem]=find(IdM==0);

idFreeParticle=[idFreeParticle;idFreeParticleTem];
idToFollow=[idToFollow;idToFollowTem];
end

nPair=numel(idFreeParticle);

Vpair=zeros(nPair,1);

for iPair=1:nPair
   Vpair(iPair)=sum((mod.mod{i_mod(1)}.var.coord(idFreeParticle(iPair),:)-mod.mod{i_mod(2)}.var.coord(idToFollow(iPair),:)).^2,2);
end

[Cunq,~,iCunq]=unique(idFreeParticle,'row');

Nunq=size(Cunq,1);

Vunq=inf(Nunq,1);

idMin=zeros(Nunq,1);

for iPair=1:nPair
    if Vunq(iCunq(iPair)) > Vpair(iPair)
        Vunq(iCunq(iPair))=Vpair(iPair);
        idMin(iCunq(iPair))=iPair;
    end
end
%%
for iUnq=1:Nunq
    iPair=idMin(iUnq);
    p.var.coord(idFreeParticle(iPair),:)=mod.mod{i_mod(2)}.var.coord(idToFollow(iPair),:);
end
idConnect=[idFreeParticle(idMin),idToFollow(idMin)];