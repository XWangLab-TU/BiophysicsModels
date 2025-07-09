function [M,CutOff] = TimeEval(dyn,M,Fname,varargin)
%--------------------------------------------------------------------------
        % TimeEval computes the time evolution of a given @Model M in a
        % general manner. Currently Langevin and rotational Langevin
        % Equations are supported.
        % input: 
        % dyn - @dynamics object, defines all required options and parameters
        % M - @Model object, all the objects in M are considered
        % Fname - name of forces needed to alter the objects
        % optional:
        % see variable arguments
        %   See also rotation 
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dyn', @(x) isa(x,'dynamics'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('Fname', @(x) iscell(x));
ip.addParameter('tCutoff', inf, @isnumeric); %specific time duration to cut off the evolution
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('sMax', 0.01, @isnumeric); %maximal distance to move coordinates in M
ip.addParameter('iPar', [], @isnumeric);
ip.parse(dyn,M,Fname,varargin{:});
%--------------------------------------------------------------------------
nF=numel(Fname);
%--------------------------------------------------------------------------
nS=dyn.pm.n_step;
tCutoff=ip.Results.tCutoff;
sMax=ip.Results.sMax;
%--------------------------------------------------------------------------
%return a structure to indicate the senario of cutoff, nS: total steps
%done; forcedOff: due to defects, not defined cutoff conditions
CutOff=struct('nS',nS,'forcedOff',false,'dataStored',[]);
%--------------------------------------------------------------------------
%%
ftot=cell(M.n_mod,1); %total force on each object in M
Vrec=cell2struct(cell(nF,1),Fname); %potential of each force recorded 
drBYdt=cell(M.n_mod,1); %dr/dt for temporarily evaluate appropriate adaptive time step
daBYdtRot=cell(M.n_mod,1); %da/dt (in rotation) for temporarily evaluate appropriate adaptive time step
drBYdtRot=cell(M.n_mod,1); %dr/dt (in rotation) for temporarily evaluate appropriate adaptive time step
for iF=1:nF
    Vrec.(Fname{iF})=zeros(nS,1);
end
t=0;
breakOff=false;
if ~isempty(dyn.field)
    fieldOn=true;
else
    fieldOn=false;
end
for iS=1:nS
%     fprintf('iS %d  %d\n',iS,ip.Results.iPar);
%     disp([max(max(M.mod{M.i_mod.ModMembrane}.var.j_T)), M.mod{M.i_mod.ModMembrane}.var.n_coord]);
    for iF=1:nF
        if dyn.iFon(iF)==true
            updateForceOnly=true;
            imForceMod=nan;
            for im=1:M.n_mod
            if (dyn.updateModForce(im)==true) && (strcmp(M.name{im},Fname{iF})==true)
                updateForceOnly=false;
                imForceMod=im;
                break;
            end
            end
            if updateForceOnly==true
               [M.TypForce] = M.TypForce.(Fname{iF})(M); %compute each required force
            else
               [M.TypForce,M.mod{M.i_mod.(M.name{imForceMod})}] = M.TypForce.(Fname{iF})(M); %compute single modular force and update the module  
            end
        Vrec.(Fname{iF})(iS)=M.TypForce.int_V.(Fname{iF}); %record each potential energy
        end
    end
    %----------------------------------------------------------------------
    % field forces
    if fieldOn==true
    for im=1:M.n_mod
        if ~isempty(dyn.field{im})
        [M.TypForce,~] = ModField(M.TypForce,M,dyn.field{im});
        end
    end
    end
    %----------------------------------------------------------------------
    %match and add each force to corresponding objects, for example,
    %M.Typforce.ModMembrane_ModSubstrate will be added to both the forces on
    %ModMembrane and ModSubstrate, defined by index {im}
    for im=1:M.n_mod
        if ~isempty(dyn.matchFtoMod{im})
        if dyn.dynTyp(im)==0 %Langevin
            ftot{im}=zeros(M.mod{M.i_mod.(M.name{im})}.var.n_coord,3);
            nMatch=size(dyn.matchFtoMod{im},1);
            for iMatch=1:nMatch
                iF=dyn.matchFtoMod{im}(iMatch,2);
                ftot{im}=ftot{im}+M.TypForce.int_tot.(Fname{iF}){dyn.matchFtoMod{im}(iMatch,1)};
            end
        elseif dyn.dynTyp(im)==1 %rotation
            ftot{im}=zeros(M.mod{M.i_mod.(M.name{im})}.var.n_coord,6);
            nMatch=size(dyn.matchFtoMod{im},1);
            for iMatch=1:nMatch
                iF=dyn.matchFtoMod{im}(iMatch,2);
                ftot{im}=ftot{im}+M.TypForce.int_tot.(Fname{iF}){dyn.matchFtoMod{im}(iMatch,1)};
            end
        end
        end
        if fieldOn==true
        if ~isempty(dyn.field{im})
           ftot{im}=ftot{im}+M.TypForce.int_field{dyn.field{im}.i_mod};
        end
        end
    end
    %----------------------------------------------------------------------
    for im=1:M.n_mod
        if dyn.needFollow(im)==true
            nameV=M.mod{M.i_mod.(M.name{im})}.var.follow.nameV;
            nameIdNeighbor=M.mod{M.i_mod.(M.name{im})}.var.follow.nameIdNeighbor;
            idModFollower=M.mod{M.i_mod.(M.name{im})}.var.follow.idModFollower;
            idModFollowee=M.mod{M.i_mod.(M.name{im})}.var.follow.idModFollowee;
%             if iS==1
%                 updateTem=false;
%             elseif breakOff==true
%                 updateTem=false;
%             else
%                 updateTem=true;
%             end
            updateTem=false;
            [M,ftot] = dyn.follow(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,'update',updateTem);
        end
    end
    %----------------------------------------------------------------------
    dt=dyn.pm.dt_const;
    for im=1:M.n_mod
        if ~isempty(dyn.matchFtoMod{im})
        if (dyn.ifVarDt(im)==true) && (dyn.dynTyp(im)>=0)
            [varDt] = dyn.(['varDt_' M.name{im}])(M.TypForce,ftot{im},M.mod{M.i_mod.(M.name{im})});
            if varDt<dt
                dt=varDt;
            end
        end
        end
    end
    %----------------------------------------------------------------------
    for im=1:M.n_mod
        if ~isempty(dyn.matchFtoMod{im})
        if dyn.dynTyp(im)==0 %Langevin-------------------------------------
            drBYdt{im}=ftot{im}*M.mod{M.i_mod.(M.name{im})}.pm.mu;
            varDt=min(abs(sMax./vecnorm(drBYdt{im},2,2)));
            if varDt<dt
                dt=varDt;
            end
        elseif dyn.dynTyp(im)==1 %rotation---------------------------------
            [varDt,daBYdtRot{im},drBYdtRot{im}] = dynamics.rotationDt(M.mod{M.i_mod.(M.name{im})},ftot{im},dt,0,'sMax',sMax);
            if varDt<dt
                dt=varDt;
            end
        end
        end 
    end
    for im=1:M.n_mod
        if ~isempty(dyn.matchFtoMod{im})
        if dyn.dynTyp(im)==0 %Langevin-------------------------------------
           M.mod{M.i_mod.(M.name{im})}.var.coord=...
               M.mod{M.i_mod.(M.name{im})}.var.coord+drBYdt{im}*dt;
        elseif dyn.dynTyp(im)==1 %rotation---------------------------------
            [M.mod{M.i_mod.(M.name{im})}] = dynamics.rotation(M.mod{M.i_mod.(M.name{im})},ftot{im},dt,0,'sMax',sMax,...
                                                                'da',daBYdtRot{im},'dr',drBYdtRot{im});
        end
        end
    end
    %----------------------------------------------------------------------
    breakOff=false;
    for im=1:M.n_mod
        if dyn.needBreakOff(im)==true
%             edgSave=M.mod{M.i_mod.(M.name{im})}.var.edge_all;
%             coordSave=M.mod{M.i_mod.(M.name{im})}.var.coord;
%             [~,m,~] = ModMembrane(M.TypForce,M);
            [M,breakOffInfo] = dyn.(['breakOff_' M.name{im}])(M);
            if breakOffInfo.breakOff==true
                breakOff=true;
%                 disp(min(M.mod{M.i_mod.ModMembrane}.var.CanSplitMerge));pause;
            end
%             if breakOffInfo.breakOff==true
%                 nTem=numel(breakOffInfo.id);
%                 kHtem=[m.var.f.kH(edgSave(breakOffInfo.id,1)),m.var.f.kH(edgSave(breakOffInfo.id,2))];
%                 coordTem=[coordSave(edgSave(breakOffInfo.id,1),:),coordSave(edgSave(breakOffInfo.id,2),:)];
%                 CutOff.dataStored=[CutOff.dataStored;[iS*ones(nTem,1),breakOffInfo.id,kHtem,coordTem]];
%             end
        end  
    end
    %----------------------------------------------------------------------
    t=t+dt;
    if t>tCutoff
        CutOff.nS=iS;
        break;
    end
    %----------------------------------------------------------------------
    for im=1:M.n_mod
       if ~isempty(M.mod{M.i_mod.(M.name{im})}.failInfo)
           CutOff.forcedOff=true;
           break;
       end
    end
    %----------------------------------------------------------------------
    %%
%     hold on;
%     im=1;
%     quiver3(M.mod{M.i_mod.(M.name{im})}.var.coord(:,1),...
%             M.mod{M.i_mod.(M.name{im})}.var.coord(:,2),...
%             M.mod{M.i_mod.(M.name{im})}.var.coord(:,3),...
%             ftot{im}(:,1),...
%             ftot{im}(:,2),...
%             ftot{im}(:,3),'linewidth',2);
end
if ip.Results.plot_or_not==true
    figure();
    subplot(1,2,1);
    VrecTot=0;
    for iF=1:nF
       VrecTot=VrecTot+Vrec.(Fname{iF});
    end
    plot(VrecTot);
    subplot(1,2,2);
    for iF=1:nF
       plot(zscore(Vrec.(Fname{iF}))); hold on;
    end
    legend(fieldnames(Vrec));
end