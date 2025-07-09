function [M,ftot] = follow(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('M', @(x) isa(x,'model'));
            ip.addRequired('nameV', @(x) ischar(x));
            ip.addRequired('nameIdNeighbor', @(x) ischar(x));
            ip.addRequired('idModFollower', @(x) isnumeric(x));
            ip.addRequired('idModFollowee', @(x) isnumeric(x));
            ip.addRequired('ftot', @(x) iscell(x));
            ip.addParameter('update', false, @islogical);
            ip.addParameter('sMax', 0.01, @isnumeric);
            ip.parse(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,varargin{:});
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
idCoordFollower=(1:M.mod{idModFollower}.var.n_coord)';
Ner=numel(idCoordFollower);
update=ip.Results.update;
if update==false
    idCoordFollowee=zeros(M.mod{idModFollower}.var.n_coord,1);
else
    idCoordFollowee=M.mod{idModFollower}.var.follow.idCoordFollowee;
end

followExtra=false;
if isa(M.mod{idModFollowee},'ModMembrane')
    idExtra=M.mod{idModFollowee}.var.j_T;
    followExtra=true;
end
% followExtra=false; %XXX

if update==true
    M.mod{idModFollower}.var.coord=M.mod{idModFollowee}.var.coord(idCoordFollowee,:);
    for i=1:M.mod{idModFollower}.var.n_coord
    ftot{idModFollowee}(idCoordFollowee(i),:)=ftot{idModFollowee}(idCoordFollowee(i),:)+ftot{idModFollower}(i,:);
    if followExtra==true
        idExtraAtVer=idExtra(idCoordFollowee(i),:);
        idExtraAtVer=idExtraAtVer(~isnan(idExtraAtVer));
        ftot{idModFollowee}(idExtraAtVer,:)=ftot{idModFollowee}(idExtraAtVer,:)+ftot{idModFollower}(i,:);
    end
    end
else
%see if moving the follower to the neighbors of the followee would lower
%the given potential
for i=1:Ner
%     disp(i);
    d=sum((M.mod{idModFollower}.var.coord(idCoordFollower(i),:)-M.mod{idModFollowee}.var.coord).^2,2);
    % here, only allow curvy vertices to be followed, avoid vertices at the
    % ceiling connecting to clathrin
    idTem=(M.mod{idModFollowee}.var.coord(:,3)>2) & (M.mod{idModFollowee}.var.f.kH<0.2);
    d(idTem)=inf;
    
    [~,M.mod{idModFollower}.var.follow.idCoordFollowee(i)]=min(d);
    M.mod{idModFollower}.var.coord(i,:)=M.mod{idModFollowee}.var.coord(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
    ftot{idModFollowee}(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:)=ftot{idModFollowee}(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:)...
                                                                             +ftot{idModFollower}(i,:);
    if followExtra==true
        idExtraAtVer=idExtra(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
        idExtraAtVer=idExtraAtVer(~isnan(idExtraAtVer));
        ftot{idModFollowee}(idExtraAtVer,:)=ftot{idModFollowee}(idExtraAtVer,:)+ftot{idModFollower}(i,:);
    end
end
end
%--------------------------------------------------------------------------
% only ModClathrin_ModFreeParticle supported
idC_idFP=M.TypForce.otherInfo.ModClathrin_ModFreeParticle.idC_idFP;
nPair=size(idC_idFP,1);
Dpair=zeros(nPair,1);
coordFoot3=getFoot(M.mod{M.i_mod.ModClathrin});
coordFP=M.mod{M.i_mod.ModFreeParticle}.var.coord;

for iPair=1:nPair
   Dpair(iPair)=sum((coordFoot3(idC_idFP(iPair,2),:,idC_idFP(iPair,1))-coordFP(idC_idFP(iPair,3),:)).^2,2);
end
for iPair=1:nPair
    iFP=idC_idFP(iPair,3);
    IdNeighbor=M.mod{idModFollowee}.var.(nameIdNeighbor)(M.mod{idModFollower}.var.follow.idCoordFollowee(iFP),:);
    IdNeighbor=IdNeighbor(~isnan(IdNeighbor));
    nNeighbor=numel(IdNeighbor);
    coordTry=M.mod{idModFollowee}.var.coord(IdNeighbor,:);
    D_alt=zeros(nNeighbor,1);
    for j=1:nNeighbor     
       D_alt(j)=sum((coordFoot3(idC_idFP(iPair,2),:,idC_idFP(iPair,1))-coordTry(j,:)).^2,2);
    end
    [Dmin,idMin]=min(D_alt);
    if (Dmin<Dpair(iPair))
        M.mod{idModFollower}.var.coord(idCoordFollower(iFP),:)=coordTry(idMin,:);
        M.mod{idModFollower}.var.follow.idCoordFollowee(i)=IdNeighbor(idMin);
    end
end
%%
% scatter3(coordTry(:,1),coordTry(:,2),coordTry(:,3),20,'filled','MarkerFaceColor',[1 0 1]);
% R=M.mod{idModFollower}.var.coord(idCoordFollower(iFP),:);
% hold on;scatter3(R(:,1),R(:,2),R(:,3),20,'filled','MarkerFaceColor',[0 1 0]);

% %old:
% [~,V_tot]=M.TypForce.(nameV)(M);
% for i=1:Ner
%     IdNeighbor=M.mod{idModFollowee}.var.(nameIdNeighbor)(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
%     nNeighbor=numel(IdNeighbor);
%     V_alt=zeros(nNeighbor,1);
%     coordTry=zeros(nNeighbor,3);
%     coordSave=M.mod{idModFollower}.var.coord(idCoordFollower(i),:);
%     for j=1:nNeighbor
%         if ~isnan(IdNeighbor(j))
%             coordTry(j,:)=M.mod{idModFollowee}.var.coord(IdNeighbor(j),:);
%             M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=M.mod{idModFollowee}.var.coord(IdNeighbor(j),:);
%             [~,V_alt(j)]=M.TypForce.(nameV)(M,'Vonly',true);            
%             M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=coordSave;
%         else
%             V_alt(j)=inf;
%         end
%     end
%     [Vmin,idMin]=min(V_alt);
%     if (Vmin<V_tot)&&(Vmin~=0)
%         M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=coordTry(idMin,:);
%         M.mod{idModFollower}.var.follow.idCoordFollowee(i)=IdNeighbor(idMin);
%     end
% end
%--------------------------------------------------------------------------                             
end