% function [f,V_tot,identifier,idConnect] = ModClathrin_ModFreeParticle(f,mod,varargin)
function [f,V_tot] = ModClathrin_ModFreeParticle(f,mod,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('dr', 0.000001, @isnumeric);
ip.addParameter('idSub1', [], @isnumeric);
ip.addParameter('idSub2', [], @isnumeric);
ip.addParameter('Vonly', false, @islogical);
ip.addParameter('update', true, @islogical);
ip.addParameter('addOtherInfo', false, @islogical);
ip.addParameter('initOtherInfo', false, @islogical);
ip.parse(f,mod,varargin{:}); 
%----------------------------------------------------------------------------------------
dr=ip.Results.dr;
Vonly=ip.Results.Vonly;
update=ip.Results.update;
if isempty(f.otherInfo.ModClathrin_ModFreeParticle)
    update=false;
elseif isempty(f.otherInfo.ModClathrin_ModFreeParticle.idC_idFP)
    update=false;
end
addOtherInfo=ip.Results.addOtherInfo;
initOtherInfo=ip.Results.initOtherInfo;
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModFreeParticle];  %1: clathrin, 2: AP2 as free particle
%----------------------------------------------------------------------------------------
idSub1=ip.Results.idSub1;
idSub2=ip.Results.idSub2;
%----------------------------------------------------------------------------------------
%%
if isempty(idSub1)
    idSub1=(1:mod.mod{mod.i_mod.ModClathrin}.var.n_coord)';
end
if isempty(idSub2)
    idSub2=(1:mod.mod{mod.i_mod.ModFreeParticle}.var.n_coord)';
end
nClathrin=numel(idSub1);
nCtot=nClathrin*3; %3 legs per Clathrin
coordFoot3=getFoot(mod.mod{mod.i_mod.ModClathrin});
coordFoot3sub=coordFoot3(:,:,idSub1);
coordFoot=[];
for i=1:nClathrin
    coordFoot=cat(1,coordFoot,coordFoot3sub(:,:,i));
end
nFP=numel(idSub2);
coordFP=mod.mod{mod.i_mod.ModFreeParticle}.var.coord;
coordFPsub=coordFP(idSub2,:);
%----------------------------------------------------------------------------------------
if update==false
%----------------------------------------------------------------------------------------
idClathrin=zeros(nCtot,2);
idClathrin(1:3:nCtot,1)=idSub1;
idClathrin(2:3:nCtot,1)=idSub1;
idClathrin(3:3:nCtot,1)=idSub1;
idClathrin(1:3:nCtot,2)=1;
idClathrin(2:3:nCtot,2)=2;
idClathrin(3:3:nCtot,2)=3;
%----------------------------------------------------------------------------------------
identifier=zeros(mod.n_mod,1);
identifier(mod.i_mod.ModClathrin)=1;
identifier(mod.i_mod.ModFreeParticle)=2;
%----------------------------------------------------------------------------------------
idC_idFP=[];
for i=1:nCtot
    D=sum((coordFoot(i,:)-coordFPsub).^2,2);
    [~,idMin]=min(D);
    idC_idFP=[idC_idFP;idClathrin(i,1),idClathrin(i,2),idSub2(idMin)];
end
otherInfo=struct('idC_idFP',idC_idFP,'identifier',identifier);
if initOtherInfo==true
    f.otherInfo.ModClathrin_ModFreeParticle=struct('idC_idFP',[],'identifier',identifier);
end
if addOtherInfo==true
f.otherInfo.ModClathrin_ModFreeParticle.idC_idFP=[f.otherInfo.ModClathrin_ModFreeParticle.idC_idFP;otherInfo.idC_idFP];
end
%----------------------------------------------------------------------------------------
else
idC_idFP=f.otherInfo.ModClathrin_ModFreeParticle.idC_idFP;
end

f.int_comp.ModClathrin_ModFreeParticle=cell(2,1);
f.int_comp.ModClathrin_ModFreeParticle{1}=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
f.int_comp.ModClathrin_ModFreeParticle{2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);

nPair=size(idC_idFP,1);
Vpair=zeros(nPair,1);
Fpair=zeros(nPair,3);

for iPair=1:nPair
   Vpair(iPair)=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((coordFoot3(idC_idFP(iPair,2),:,idC_idFP(iPair,1))-coordFP(idC_idFP(iPair,3),:)).^2,2);
   Fpair(iPair,:) = -f.pm.k_ModClathrin_ModFreeParticle*(coordFP(idC_idFP(iPair,3),:)-...
                                               coordFoot3(idC_idFP(iPair,2),:,idC_idFP(iPair,1)));
end

% idConnect=[idClathrin(idMin,:),idFreeParticle(idMin)];
V_tot=sum(Vpair);
f.int_V.ModClathrin_ModFreeParticle=V_tot;
%%
if Vonly==false
id_feet=[49 50 51];
%%
% tic
for iPair=1:nPair
    for m=1:3
%--------------------------------------------------------------------------
       coordFootTem=coordFoot3; 
       Vtem=V_tot-Vpair(iPair);
       coordFootTem(idC_idFP(iPair,2),m,idC_idFP(iPair,1))=coordFootTem(idC_idFP(iPair,2),m,idC_idFP(iPair,1))+dr;
%--------------------------------------------------------------------------   
       VpairTem=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((coordFootTem(idC_idFP(iPair,2),:,idC_idFP(iPair,1))-coordFP(idC_idFP(iPair,3),:)).^2,2);
       Vtem=Vtem+VpairTem;
       f.int_comp.ModClathrin_ModFreeParticle{1}(idC_idFP(iPair,1),m)=...
       f.int_comp.ModClathrin_ModFreeParticle{1}(idC_idFP(iPair,1),m)-(Vtem-V_tot)/dr;
       f.int_comp.ModClathrin_ModFreeParticle{2}(idC_idFP(iPair,3),m)=...
       f.int_comp.ModClathrin_ModFreeParticle{2}(idC_idFP(iPair,3),m)+(Vtem-V_tot)/dr;
%--------------------------------------------------------------------------
       c=mod.mod{i_mod(1)};
       Vtem=V_tot-Vpair(iPair);
%--------------------------------------------------------------------------       
       c.var.a(idC_idFP(iPair,1),m)=c.var.a(idC_idFP(iPair,1),m)+dr;
       c.var.ang_a.Phi(idC_idFP(iPair,1))=norm(c.var.a(idC_idFP(iPair,1),:));
%        [c.var] = setVar(c,true); %below 4 lines work the same but faster
       var=c.var;
       var.O=c.Omega(var.a,var.ang_a.Phi);
       var.E=c.getE(var.a,var.ang_a.Phi);
       c.var=var;
       
       foot_tem=c.get_r_from_a(c.var.coord_org(id_feet(idC_idFP(iPair,2)),1:3),1,c.var.O(:,:,idC_idFP(iPair,1)))+c.var.coord(idC_idFP(iPair,1),:);    
       VpairTem=0.5*f.pm.k_ModClathrin_ModFreeParticle*...
          sum((foot_tem-coordFP(idC_idFP(iPair,3),:)).^2,2);
       Vtem=Vtem+VpairTem;
       f.int_comp.ModClathrin_ModFreeParticle{1}(idC_idFP(iPair,1),m+3)=...
              f.int_comp.ModClathrin_ModFreeParticle{1}(idC_idFP(iPair,1),m+3)-(Vtem-V_tot)/dr;
%--------------------------------------------------------------------------
    end
end
% toc
f.int_tot.ModClathrin_ModFreeParticle=f.int_comp.ModClathrin_ModFreeParticle;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end


