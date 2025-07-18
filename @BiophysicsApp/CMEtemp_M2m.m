dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
cd(dirRoot);
%% ========================================================================
%% ========================================================================
nCLmax=25;
nParaX=2;
nParaY=10;
nPara=nParaX*nParaY;
nTdef=10000;
nRep=1;
nPar=12;
xyzLimPlot=[-15 15;-5 25;-26.2 3.8];viewAng=[0 0];
rMaxFromCtr=10;
rMinFromCtr=5;
Tcutoff=5.0000e-08;
Machine='M2';
initMem=[3 5];

dirFolder=cell(nCLmax,1);
temFormat=['%0.' num2str(3) 'd'];
for iCL=1:nCLmax
    dirFolder{iCL}=[dirRoot filesep Machine filesep num2str(iCL,temFormat)]; mkdir(dirFolder{iCL});
    rec=ComRecord([dirFolder{iCL} filesep 'data'],[dirFolder{iCL} filesep 'fig'],[dirFolder{iCL} filesep 'init'],...
              'nRep',nRep,'deleteFiles',true,...
              'paraAlt',struct('kappa',[30 50],'kCL',(20:5:65))); %(20:10:90)
    save([dirFolder{iCL} filesep 'rec.mat'],'rec');
end

nTeq=nTdef*ones(nCLmax,nPara);
nTeq(1,:)=-1;

secDone=false(nCLmax,nPara);
secDone(1,:)=true;

%%
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_new_test_withCL.mat');
M1=S.Mpar{initMem(1)};
M2=S.Mpar{initMem(2)};
clear S;
S=load([dirFolder{1} filesep 'rec.mat']);
i_loop_done=0;
for iY=1:nParaY
for iX=1:nParaX
    iFolder=(iY-1)*nParaX+iX;
    if iX==1
        M=M1;
    elseif iX==2
        M=M2;
    end
    M.TypForce.pm.k_ModClathrin_ModFreeParticle=100;
    M.TypForce.pm.k_ModClathrin(1)=rec.dir_all{iFolder}.paraVal(2); %60
    M.TypForce.pm.k_ModClathrin(2)=200.; %1000
    M.mod{M.i_mod.ModClathrin}.pm.n_min_adp=3;
    M.mod{M.i_mod.ModClathrin}.pm.n_max_adp=3;
    M.mod{M.i_mod.ModFreeParticle}.prop = {'Particle','needFollow'};
%     M.mod{M.i_mod.ModMembrane}.pm.DlocRelax=20;
%     M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.5+0.*iPar;
%     M.mod{M.i_mod.ModMembrane}.pm.mu=500000; %200000
%     M.mod{M.i_mod.ModMembrane}.pm.k_c=100;
%     M.TypForce.int_stored.ModMembrane=[];
    dirTem=S.rec.dir_all{iFolder};
    fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirTem,M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
end
end
clear S;
%%
for iCL=2:nCLmax
    S=load([dirFolder{iCL} filesep 'rec.mat']);
for iPara=1:nPara
    dirTem=S.rec.dir_all{iPara};
    iFrames = ComRecord.getFrame(dirTem.dir_data);
    if isempty(iFrames)
        iFrames=0;
    end
    if max(iFrames) >= nTeq(iCL,iPara)
        secDone(iCL,iPara)=true;
    else
        secDone(iCL,iPara)=false;
    end
end
    clear S;
end
%
toDo=[];
for iCL=2:nCLmax
    for iPara=1:nPara
        if (secDone(iCL,iPara)==false) && (secDone(iCL-1,iPara)==true)
            toDo=[toDo; [iCL,iPara]];
        end
    end
end
NtoDo=size(toDo,1);
if NtoDo>nPar
    nParToDo=nPar;
    ItoDo=randsample(NtoDo,nPar);
else
    nParToDo=NtoDo;
    ItoDo=(1:NtoDo)';
end
%
dirFolderPar=cell(nPar,1);
addCLorNot=cell(nPar,1);
Mpar=cell(nPar,1);
i_loop_done_all=cell(nPar,1);
nTeqPar=cell(nPar,1);
iClathrin=cell(nPar,1);
CtrPar=cell(nPar,1);
for iPar=1:nParToDo
    iCL=toDo(ItoDo(iPar),1);
    iPara=toDo(ItoDo(iPar),2);
    S=load([dirFolder{iCL} filesep 'rec.mat']);
    dirFolderPar{iPar}=S.rec.dir_all{iPara};
    dirTem=dir(dirFolderPar{iPar}.dir_data);
    if numel(dirTem) <=2
        addCLorNot{iPar}=true;  
        s=load([dirFolder{iCL-1} filesep 'rec.mat']);
        dirTem=dir(s.rec.dir_all{iPara}.dir_data);
        ss=load([dirTem(end).folder filesep dirTem(end).name]);
        Mpar{iPar}=ss.M;
        i_loop_done_all{iPar}=0;
        clear s;
    else
        addCLorNot{iPar}=false; 
        dirTem=dir(S.rec.dir_all{iPara}.dir_data);
        ss=load([dirTem(end).folder filesep dirTem(end).name]);
        Mpar{iPar}=ss.M;
        iFrames=ComRecord.getFrame(S.rec.dir_all{iPara}.dir_data);
        i_loop_done_all{iPar}=max(iFrames);
    end
    nTeqPar{iPar}=nTeq(iCL,iPara);
    iClathrin{iPar}=iCL;
    CtrPar{iPar}=[0 0 0];
end
clear S;
clear ss;
% add clathrin and relax
save([dirRoot filesep Machine filesep 'progress0.mat']);
parfor iPar=1:nParToDo
% for iPar=4:4
    %%---------------------------------------------------------------------
    M=Mpar{iPar};
    M.mod{M.i_mod.ModMembrane}.pm.n_val_min=5;
    M.mod{M.i_mod.ModMembrane}.pm.n_val_max=8;
%     M.mod{M.i_mod.ModMembrane}.var.CanSplitMerge=ones(M.mod{M.i_mod.ModMembrane}.var.n_edg,2);
%     M.mod{M.i_mod.ModMembrane}.pm.remeshExtraCtrl=1000;
%     M.mod{M.i_mod.ModMembrane}.pm.doRemeshExtraCtrl=true;
    M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=1;
    M.mod{M.i_mod.ModMembrane}.pm.DlocRelax=100;
    %%---------------------------------------------------------------------
    if addCLorNot{iPar}==true
        if iClathrin{iPar}<4
            near=true;
        else
            near=false;
        end
        fprintf('add Clathrin...\n');
        if near==true
            name={'ModClathrin_ModFreeParticle'};
            [M,~] = addNear(M.mod{M.i_mod.ModClathrin},M,CtrPar{iPar},'rMaxFromCtr',rMaxFromCtr,'xyzLim',[-inf,inf;-inf,inf;-4,0]);
            idNew1=M.mod{M.i_mod.ModClathrin}.var.n_coord;
            idNew2=(M.mod{M.i_mod.ModFreeParticle}.var.n_coord-2:M.mod{M.i_mod.ModFreeParticle}.var.n_coord);
            [M] = locDyn(M.mod{M.i_mod.ModClathrin},M,name,idNew1,idNew2,'sMax',1,'update',false,'initOtherInfo',false,'addOtherInfo',true);
        else
            %add clathrin from far distance and locate it back to the coat
            muSave=M.mod{M.i_mod.ModMembrane}.pm.mu;
            M.mod{M.i_mod.ModMembrane}.pm.mu=0.01; %temporarily immobilize the membrane
            %----------------------------------------------------------------------
            [M,~] = add(M.mod{M.i_mod.ModClathrin},M,[0 0 -3.5],'rMaxFromCtr',rMaxFromCtr,'rMinFromCtr',rMinFromCtr);
            nStepTem=1000;
            field_para=cell(1,1);
            field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',-3.25,'z_upper',4.25,'order',1,'steep',5));
            dyn=dynamics(nStepTem);
            Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
            [dyn] = Preparation(dyn,M,Fname,'field_or_not',true,'field_para',field_para);
            dyn.needBreakOff(1)=false;
            i_loop_done=0;
            while (i_loop_done<nStepTem)
                [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*0.01,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
                %dyn.breakOff_ModClathrin(M)
                %link(M.mod{M.i_mod.ModClathrin})
                i_loop_done=i_loop_done+CutOff.nS;
                fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
            end
            %----------------------------------------------------------------------
            M.mod{M.i_mod.ModMembrane}.pm.mu=muSave;
        end
    end
    %%---------------------------------------------------------------------
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    field_para=cell(1,1);
    field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',-3.25,'z_upper',4.25,'order',1,'steep',5));
    dyn=dynamics(nTeqPar{iPar});
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname,'field_or_not',true,'field_para',field_para);
    i_loop_done=i_loop_done_all{iPar};
    fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirFolderPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
    while (i_loop_done<nTeqPar{iPar})
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*0.2,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
%         f=M.TypForce;
%         eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
        i_loop_done=i_loop_done+CutOff.nS;
        fig=figure;
        ComRecord.addFrame(dirFolderPar{iPar},M,fig,i_loop_done,'rec_mod_only',true,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirFolderPar{iPar},M,fig,i_loop_done,'rec_fig_only',true,'update_mod_only',false); 
    Mpar{iPar}=M;
    i_loop_done_all{iPar}=i_loop_done;   
end

% determine whether equilibrium is reached
disp('reading results...');
PF=true(nParToDo,1);
Yplot=cell(nParToDo,1);
parfor iPar=1:nParToDo
    dirTem=dir(dirFolderPar{iPar}.dir_data);
    nFile=numel(dirTem);
    Yplot{iPar}=zeros(nFile,1);
    for iFile=4:nFile
        S=load([dirTem(iFile).folder filesep dirTem(iFile).name]);
        M=S.M;
        f=M.TypForce;
        S.M.TypForce.pm.k_ModClathrin=[1 0];
        [f,~,~] = ModClathrin(f,S.M);
        nLeg=0;
        for i=1:S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord
            nLeg=nLeg+sum(sum(S.M.mod{S.M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
        end
        nLeg=nLeg*0.5;
        x=S.M.mod{S.M.i_mod.ModClathrin}.var.coord_org;
        lPL=norm(x(52,1:3)-x(8,1:3));
        Yplot{iPar}(iFile)=1./(sqrt(f.int_V.ModClathrin)/nLeg/lPL);
    end
end
save([dirRoot filesep Machine filesep 'progress1.mat']);
for iPar=1:nParToDo
    fig=figure;
    ax = axes(fig);
    plot(ax,Yplot{iPar}(4:end),'-*');
    PF(iPar)=ComPlot.PassFail(fig,dirFolderPar{iPar}.dir_data);
end
for iPar=1:nParToDo
    iCL=toDo(ItoDo(iPar),1);
    iPara=toDo(ItoDo(iPar),2);
    if PF(iPar)==false
        nTeq(iCL,iPara)=nTeq(iCL,iPara)+nTdef;
    end
end
save([dirRoot filesep Machine filesep 'progress2.mat']);
disp('done!');
%%
dispAllFig=true;
if dispAllFig==true
    temFormat2=['%0.' num2str(2) 'd'];
    for iPara=1:nPara
        for iCL=2:nCLmax
            if ((secDone(iCL,iPara)==false) && (secDone(iCL-1,iPara)==true))
                dirTem=dir([dirFolder{iCL} filesep 'data' filesep 'p' num2str(iPara, temFormat2) '_rep1']);
                SS=load([dirTem(end).folder filesep dirTem(end).name]);
                M=SS.M;
                c=M.mod{M.i_mod.ModClathrin};
                [Ttem] = getTightness(c,M);
%                 SSS=load([dirFolder{iCL} filesep 'rec.mat']); dirPlot=SSS.rec.dir_all{iPara}.dir_fig; dirTem=dir(dirPlot);
%                 fig=openfig([dirTem(end).folder filesep dirTem(end-1).name]); 
                fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng); 
                iFrames = ComRecord.getFrame([dirTem(end).folder]);
                if isempty(iFrames)
                    iFrames=0;
                else
                    iFrames=max(iFrames);
                end
                title([num2str(iPara),'; ', num2str(iCL), '; ',num2str(Ttem), '; ' num2str(iFrames), '; ', num2str(nTeq(iCL,iPara))]);
                pause; close(fig);
            end
        end
    end
    clear SS; %clear SSS;
end
% nTeq(:,11)=-1;