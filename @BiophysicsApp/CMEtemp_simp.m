dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
cd(dirRoot);
%%
doInit=false;
Machine='M3';
Midx=str2double(Machine(2));
xyzLimPlot=[-15 15;-5 25;-26.2 3.8];viewAng=[0 0];

addDyn=false;
DynStrong=false;
k1Def=true;
adj_K3_no_db_hf=1; %0: default, 1: double, 2: half

if DynStrong==true
    kDyn=[20 1]; % 1: direction
else
    kDyn=[10 1];
end

% pathMachine=[dirRoot filesep Machine 'simp']; mkdir(pathMachine);
% pathMachine=[dirRoot filesep Machine 'simp_dyn']; mkdir(pathMachine);
% pathMachine=[dirRoot filesep Machine 'simp_dyn_weak']; mkdir(pathMachine);
% pathMachine=[dirRoot filesep Machine 'simp_pdf']; mkdir(pathMachine);
pathMachine=[dirRoot filesep Machine 'simp_K3_db']; mkdir(pathMachine);
% pathMachine=[dirRoot filesep Machine 'simp_K3_hf']; mkdir(pathMachine);
%% init
if doInit==true
% S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/M2/020/data/p10_rep1/mod_Stg_001Fr_0130848.mat');
%S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/M1/015/data/p03_rep1/mod_Stg_001Fr_0350171.mat');
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/M1/016/data/p03_rep1/mod_Stg_001Fr_0251265.mat');
M=S.M;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);

rMaxFromCtr=10;
rMinFromCtr=5;
Tcutoff=5.0000e-08;
%add clathrin from far distance and locate it back to the coat
muSave=M.mod{M.i_mod.ModMembrane}.pm.mu;
M.mod{M.i_mod.ModMembrane}.pm.mu=0.01; %temporarily immobilize the membrane
%----------------------------------------------------------------------
for iAdd=1:4
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
    [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*0.01,'plot_or_not',false,'sMax', 0.1);
    %dyn.breakOff_ModClathrin(M)
    %link(M.mod{M.i_mod.ModClathrin})
    i_loop_done=i_loop_done+CutOff.nS;
    fprintf('loop finished %d \n', i_loop_done);
end
end
%----------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.pm.mu=muSave;
end
%%
%% ======================================================================== setting dir
nCLmax=20;
iCLstart=20;
nParaX=2;
nParaY=6;
nPara=nParaX*nParaY;
nTdef=40000;
nRep=6;
nPar=12;
xyzLimPlot=[-15 15;-5 25;-26.2 3.8];viewAng=[0 0];
%rMaxFromCtr=10;
%rMinFromCtr=5;
Tcutoff=5.0000e-08;
%kCLperc=[80 85 90 95 100 105 110 115];
%kCLperc=[75+(Midx-1)*10 80+(Midx-1)*10];
if k1Def==true
    kCLperc=[60+(3-1)*20 60+(3-1)*20]; %both Midx=3 2nd half
else
    kCLperc=[50+(Midx-1)*20 60+(Midx-1)*20];
end

dirFolder=cell(nCLmax,1);
temFormat=['%0.' num2str(3) 'd'];
for iCL=iCLstart:nCLmax
    dirFolder{iCL}=[pathMachine filesep num2str(iCL,temFormat)]; mkdir(dirFolder{iCL});
    rec=ComRecord([dirFolder{iCL} filesep 'data'],[dirFolder{iCL} filesep 'fig'],[dirFolder{iCL} filesep 'init'],...
              'nRep',nRep,'deleteFiles',true,...
              'paraAlt',struct('kCLperc',kCLperc)); %(20:10:90)
    save([dirFolder{iCL} filesep 'rec.mat'],'rec');
end

nTeq=nTdef*ones(nCLmax,nPara);
nTeq(1,:)=-1;
nTeq(:,:)=-1;
nTeq(nCLmax,:)=nTdef;

secDone=false(nCLmax,nPara);
secDone(1,:)=true;

secDone(:,:)=true;
secDone(nCLmax,:)=false;
% ======================================================================== setting parameters
S=load([pathMachine filesep 'init.mat']);
Msave=S.M;
clear S;
S=load([dirFolder{nCLmax} filesep 'rec.mat']);
i_loop_done=0;
Mpar=cell(nPar,1);
for iX=1:nParaX
for iY=1:nParaY
    M=Msave;
    iFolder=(iY-1)*nParaX+iX;
    M.TypForce.pm.k_ModClathrin(1)=rec.dir_all{iFolder}.paraVal/100*Msave.TypForce.pm.k_ModClathrin(1); %60
    if adj_K3_no_db_hf==1
        M.TypForce.pm.k_ModClathrin_ModFreeParticle=Msave.TypForce.pm.k_ModClathrin_ModFreeParticle*1.5;
    elseif adj_K3_no_db_hf==2
        M.TypForce.pm.k_ModClathrin_ModFreeParticle=Msave.TypForce.pm.k_ModClathrin_ModFreeParticle*0.5;
    elseif adj_K3_no_db_hf~=0
        error('adj_K3_no_db_hf wrong!!')
    end
    if addDyn==true
        M.TypForce.pm.k_ModFreeParticle_ModMembrane=kDyn;
    end
    dirTem=S.rec.dir_all{iFolder};
    fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirTem,M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
    Mpar{iFolder}=M;
end
end
clear S;
%% ======================================================================== loop
iCL=nCLmax;
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

toDo=[];
    for iPara=1:nPara
        if (secDone(iCL,iPara)==false) && (secDone(iCL-1,iPara)==true)
            toDo=[toDo; [iCL,iPara]];
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

dirFolderPar=cell(nPar,1);
i_loop_done_all=cell(nPar,1);
nTeqPar=cell(nPar,1);
iClathrin=cell(nPar,1);
for iPar=1:nParToDo
    iPara=toDo(ItoDo(iPar),2);
    S=load([dirFolder{iCL} filesep 'rec.mat']);
    dirFolderPar{iPar}=S.rec.dir_all{iPara};
    dirTem=dir(S.rec.dir_all{iPara}.dir_data);
    ss=load([dirTem(end).folder filesep dirTem(end).name]);
    Mpar{iPar}=ss.M;
    iFrames=ComRecord.getFrame(S.rec.dir_all{iPara}.dir_data);
    i_loop_done_all{iPar}=max(iFrames);
    nTeqPar{iPar}=nTeq(iCL,iPara);
    iClathrin{iPar}=iCL;
end
clear S;
clear ss;
save([pathMachine filesep 'progress0.mat']);
%--------------------------------------------------------------------------
parfor iPar=1:nParToDo
% for iPar=4:4
    %%---------------------------------------------------------------------
    M=Mpar{iPar};
    %%---------------------------------------------------------------------
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    field_para=cell(1,1);
    field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',-3.25,'z_upper',4.25,'order',1,'steep',5));
    dyn=dynamics(nTeqPar{iPar});
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    if addDyn==false
        Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    else
        Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModFreeParticle_ModMembrane'};
    end
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
save([pathMachine filesep 'progress1.mat']);
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
        nTeq(iCL,iPara)=nTeq(iCL,iPara)+nTdef*0.8;
    end
end
save([pathMachine filesep 'progress2.mat']);
disp('done!');
% nTeq(end,:)=nTeq(end,:)+40000;