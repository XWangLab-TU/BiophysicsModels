function [] = ClathrinMediatedEndocytosisInit(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
%         ClathrinMediatedEndocytosisInit performs the initiation for
%         ClathrinMediatedEndocytosis
%         input: 
%         dirRoot - root directory for data storage
%         dirMod - directory of @ModMembrane and other required objects
%         optional:
%         see variable arguments
%           See also tetherPull, membraneFussion, exampleLamellipodia,
%           exampleEndocytosis, redBloodCell
%         Author: Xinxin Wang, Danuser Lab
%         email: wangxinxin8627@gmail.com
%         date: 2023/08/27
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', [10000 2000], @isnumeric); 
%1: simulation steps for after adding CALM control points; 
%2: simulation steps for relaxation after each clathrin added
ip.addParameter('nRep', 4, @isnumeric); %repeat number, repeat 4 times for each initial curvature: flat, shallow and high
ip.addParameter('xyzLim', [-20 20;-20 20;-20 20], @isnumeric); %mesh range
ip.addParameter('xyzLimPlot', [-10 10;-10 10;-15 15], @isnumeric); %mesh range
ip.addParameter('viewAng', [0 0], @isnumeric); %figure view angle
ip.addParameter('Tcutoff', 0.00000005, @isnumeric); %time cutoff for when to exit and plot frames
ip.addParameter('loadRes', true, @islogical); 
ip.addParameter('nCL', 20, @isnumeric);
% ip.addParameter('idCALM', [1 5 9], @isnumeric); 
ip.parse(dirRoot,dirMod,varargin{:});
%--------------------------------------------------------------------------
nRep=ip.Results.nRep;
xyzLim=ip.Results.xyzLim;
viewAng=ip.Results.viewAng;
xyzLimPlot=ip.Results.xyzLimPlot;
Tcutoff=ip.Results.Tcutoff;
loadRes=ip.Results.loadRes;
nCL=ip.Results.nCL;
%--------------------------------------------------------------------------
dirRoot=ip.Results.dirRoot;
dirMod=ip.Results.dirMod;
% dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
% xyzLimPlot=[-15 15;-5 25;-15 15];viewAng=[0 0]; xyzLim=[-20 20;-20 20;-20 20]; z_lower=-3.25;z_upper=5.75;
% Tcutoff=0.00000005; nRep=12; loadRes=false;
%--------------------------------------------------------------------------
%% initial parameters
%%
% Initialization
%--------------------------------------------------------------------------
%setup parameters and variables
% unit: 10nm, 1kBT as natural units 
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
% membrane
m=ModMembrane(true,3,0,'unit',u);
% control points for pulling membrane as init conditions
s=ModSubstrate(0.01,13);
% clathrin
c=ModClathrin('unit',u);
% later as AP2-adaptor proteins to connect clathrin and membrane
p=ModFreeParticle(0,'unit',u);
p.prop = {'Particle','needFollow'};

%introducing various interactions to be included, for example: {m,s} means ctrl points interact with membrane
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{c,p},{m,s},{p,m}};{{}}});
% assemble the model with the objects
M = model(int_info,dirMod,u,xyzLim,m,s,c,p);
%% computation
% other related parameter setting
%--------------------------------------------------------------------------
% setup path for saving results, iCALM indicates 3 different initial
% curvature: flat, shallow and high
dirFolder=[dirRoot filesep 'CALMinit']; mkdir(dirFolder);
if loadRes==true
    S=load([dirFolder filesep 'rec.mat']);
    rec=S.rec;
    clear S;
else
    rec=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],...
                'deleteFiles',true);  
%               'nRep',nRep,'deleteFiles',true,'paraAlt',struct('P',P,'kA',kA,'iCALM',[0 1 2]));              
    save([dirFolder filesep 'rec.mat'],'rec');
end
nPar=numel(rec.dir_all); %total number of workers, defult: 12
%--------------------------------------------------------------------------
% get control points' positions
[minZ,idTem]=min(M.mod{M.i_mod.ModMembrane}.var.coord(:,3));
R=abs(minZ);
S0=4*pi*R^2;
V0=4/3*pi*R^3;
% adjust membrane parameters, see ModMembrane for more info
%--------------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.pm.A0=S0*3; % here, introduce some tension to begin
M.mod{M.i_mod.ModMembrane}.pm.V0=V0*1.5;
M.mod{M.i_mod.ModMembrane}.pm.k_A=8;
M.mod{M.i_mod.ModMembrane}.pm.k_c=8;
M.mod{M.i_mod.ModMembrane}.pm.k_V=0;
M.mod{M.i_mod.ModMembrane}.pm.P=1;
M.mod{M.i_mod.ModMembrane}.pm.mu=500000; %200000
M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.5; %lipid reorganization
M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
M.mod{M.i_mod.ModMembrane}.pm.DlocRelax=20;
%--------------------------------------------------------------------------
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
%     if iPar<7
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=50;
%     else
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=100;
%     end
      Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=10*iPar;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.P=rec.dir_all{iPar}.paraVal(1);
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=rec.dir_all{iPar}.paraVal(2);
%     Mpar{iPar}.TypForce.pm.k_ModMembrane_ModSubstrate=10;  % for membrane-ctrl point force
%     if rec.dir_all{iPar}.paraVal(3)==0
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordFlat;
%     elseif rec.dir_all{iPar}.paraVal(3)==1
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordLow;
%     else
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordHigh;
%     end
    dirPar{iPar}=rec.dir_all{iPar};
end
%--------------------------------------------------------------------------
%% initial morphology
%--------------------------------------------------------------------------
parfor iPar=1:nPar
%     delete([dirPar{iPar}.dir_data filesep '*.*']);
%     delete([dirPar{iPar}.dir_fig filesep '*.*']);
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
end
fprintf('computing global morphology ...\n');
nStep=100000;
parfor iPar=1:nPar
% for iPar=1:1
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModMembrane'};
    field_para=cell(1,1);
    field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',z_lower,'z_upper',z_upper,'order',1,'steep',5));
    [dyn] = Preparation(dyn,M,Fname,'field_or_not',true,'field_para',field_para);
    dyn.needFollow(2)=false; %avoid computing AP2 for now
    i_loop_done=0;
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'plot_or_not',false,'sMax', 0.1);
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end 
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('inital morphology finishded at i_par %d \n',iPar);
end
for iPar=1:nPar
   dirTem=dir(rec.dir_all{iPar}.dir_data);
   S=load([dirTem(end).folder filesep dirTem(end).name]);
   Mpar{iPar}=S.M;
end
save('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_new_test_memOnly.mat');
%% add clathrin
dirFolder=[dirRoot filesep '1stCL']; mkdir(dirFolder);
nRep=12;
if loadRes==true
    S=load([dirFolder filesep 'rec.mat']);
    recCL=S.rec;
    clear S;
else
    recCL=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],...
              'nRep',nRep,'deleteFiles',true);
    save([dirFolder filesep 'rec.mat'],'rec');
end
nPar=numel(recCL.dir_all); %total number of workers, defult: 12  
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_new_test_memOnly.mat');
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    dirPar{iPar}=recCL.dir_all{iPar};
    Mpar{iPar}=S.Mpar{iPar};
end
clear S;
CtrPar=cell(nPar,1);
nCtr=1;
for iPar=1:nPar
    R=Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.var.coord;
    idDown=R(:,3)<0;
    R=R(idDown,:);
    [~,idSort]=sort(vecnorm(R(:,1:2),2,2));
    Zctr=abs(R(idSort(1),3));
    
    CtrPar{iPar}=zeros(nCtr,3);
    CtrPar{iPar}(1,:)=[0,0,-(Zctr-3)];
%     CtrPar{iPar}(2,:)=[0,0,(Zctr-3)];
end
fprintf('initializing 1st Clathrin...\n');
parfor iPar=1:nPar
% for iPar=1:1 
    M=Mpar{iPar};
    for iCtr=1:nCtr
            if iCtr==1
                update=false;
            else
                update=true;
            end
    changed=false;
    while (changed==false)
        [M,changed] = init(M.mod{M.i_mod.ModClathrin},M,CtrPar{iPar}(iCtr,:),'update',update,'xyzLim',[-inf,inf;-inf,inf;-4,0]);
    end
    end
    Mpar{iPar}=M;
end

name={'ModClathrin_ModFreeParticle'};
parfor iPar=1:nPar
% for iPar=1:1  
    M=Mpar{iPar};
    [M] = locDyn(M.mod{M.i_mod.ModClathrin},M,name,1,(1:3),'sMax',0.01,'update',false,'initOtherInfo',true,'addOtherInfo',true);
    Mpar{iPar}=M;
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
end
save('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_new_test_withCL.mat');
%%
xyzLimPlot=[-15 15;-5 25;-27.2 2.8];viewAng=[0 0];
for iPar=1:nPar
    dirTem=dir(dirPar{iPar}.dir_data);
    S=load([dirPar{iPar}.dir_data filesep dirTem(end).name]);
    Mpar{iPar}=S.M;
end
clear S;
%%
nStep=200000;
renew=true;
nPar=12;
rMaxFromCtr=10;
rMinFromCtr=6;
% renew=false; nStep=35000;
for iClathrin=21:21%nCL
    %%---------------------------------------------------------------------
    dirFolder=[dirRoot filesep 'CLadding' filesep num2str(iClathrin)]; mkdir(dirFolder);
    recAdd=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],'deleteFiles',true);    
    for iPar=1:nPar
        dirPar{iPar}=recAdd.dir_all{iPar};
    end
    %%---------------------------------------------------------------------
%     if iClathrin<=6
%         rMaxFromCtr=5; %5(6), 6(10), 6.5(13)
%     elseif iClathrin<=10
%         rMaxFromCtr=6; 
%     elseif iClathrin<=13
%         rMaxFromCtr=6.5; 
%     else
%         rMaxFromCtr=7.5; 
%     end
    if iClathrin<8
        near=true;
    else
        near=false;
    end
%%
fprintf('add Clathrin...\n');
if near==true
name={'ModClathrin_ModFreeParticle'};
parfor iPar=1:nPar
% for iPar=1:1  
    M=Mpar{iPar};
    for iCtr=1:nCtr
       [M,~] = addNear(M.mod{M.i_mod.ModClathrin},M,CtrPar{iPar}(iCtr,:),'rMaxFromCtr',rMaxFromCtr,'xyzLim',[-inf,inf;-inf,inf;-4,0]);
       idNew1=M.mod{M.i_mod.ModClathrin}.var.n_coord;
       idNew2=(M.mod{M.i_mod.ModFreeParticle}.var.n_coord-2:M.mod{M.i_mod.ModFreeParticle}.var.n_coord);
       [M] = locDyn(M.mod{M.i_mod.ModClathrin},M,name,idNew1,idNew2,'sMax',1,'update',false,'initOtherInfo',false,'addOtherInfo',true);
    end
    Mpar{iPar}=M;
end
else
%--------------------------------------------------------------------------
%add clathrin from far distance and locate it back to the coat
parfor iPar=1:nPar
    M=Mpar{iPar};
    muSave=M.mod{M.i_mod.ModMembrane}.pm.mu;
    M.mod{M.i_mod.ModMembrane}.pm.mu=0.01; %temporarily immobilize the membrane
    %----------------------------------------------------------------------
    [M,~] = add(M.mod{M.i_mod.ModClathrin},M,[0 0 -3.5],'rMaxFromCtr',rMaxFromCtr,'rMinFromCtr',rMinFromCtr);
    nStepTem=2000;
    field_para=cell(1,1);
    field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',z_lower,'z_upper',z_upper,'order',1,'steep',5));
    dyn=dynamics(nStepTem);
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname,'field_or_not',true,'field_para',field_para);
    i_loop_done=0;
    while (i_loop_done<nStepTem)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*0.01,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
        i_loop_done=i_loop_done+CutOff.nS;
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    %----------------------------------------------------------------------
    M.mod{M.i_mod.ModMembrane}.pm.mu=muSave;
    Mpar{iPar}=M;
end
%--------------------------------------------------------------------------
end
%%
fprintf('relaxing. ..\n');
if renew==true
eStore=cell(nPar,1);
i_loop_done_all=cell(nPar,1);
end
for iPar=1:nPar
    if renew==true
    eStore{iPar}=[];
    i_loop_done_all{iPar}=0;
    end
end
parfor iPar=1:nPar
% for iPar=1:1
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    field_para=cell(1,1);
    field_para{1}=struct('i_mod',M.i_mod.ModMembrane,'type','stripe_z','spec',struct('z_lower',-3.25,'z_upper',4.25,'order',1,'steep',5));
    dyn=dynamics(nStep);
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname,'field_or_not',true,'field_para',field_para);
    i_loop_done=i_loop_done_all{iPar};
%     i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*0.2,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
        f=M.TypForce;
        eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    Mpar{iPar}=M;
    i_loop_done_all{iPar}=i_loop_done;
end
save('/endosome/work/bioinformatics/s171152/data/CMEmodel/CLadding/temp.mat');
end
%%
for iPar=1:nPar
    figure;
    for i=1:3
    subplot(1,3,i); plot(eStore{iPar}(:,i)); title(num2str(iPar)); hold on;
    end
end
%% ========================================================================
%% ========================================================================
nCLmax=25;
nParaX=2;
nParaY=8;
nPara=nParaX*nParaY;
nTdef=10000;
nRep=1;
nPar=12;
xyzLimPlot=[-15 15;-5 25;-26.2 3.8];viewAng=[0 0];
rMaxFromCtr=10;
rMinFromCtr=5;
Tcutoff=5.0000e-08;

dirFolder=cell(nCLmax,1);
temFormat=['%0.' num2str(3) 'd'];
for iCL=1:nCLmax
    dirFolder{iCL}=[dirRoot filesep 'CLadding' filesep num2str(iCL,temFormat)]; mkdir(dirFolder{iCL});
    rec=ComRecord([dirFolder{iCL} filesep 'data'],[dirFolder{iCL} filesep 'fig'],[dirFolder{iCL} filesep 'init'],...
              'nRep',nRep,'deleteFiles',true,...
              'paraAlt',struct('kappa',[10 25],'kCL',(20:10:90)/2.5)); %(20:10:90)
    save([dirFolder{iCL} filesep 'rec.mat'],'rec');
end

nTeq=nTdef*ones(nCLmax,nPara);
nTeq(1,:)=-1;

secDone=false(nCLmax,nPara);
secDone(1,:)=true;

%%
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_new_test_withCL.mat');
M1=S.Mpar{1};
M2=S.Mpar{7};
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
    if max(iFrames) > nTeq(iCL,iPara)
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

parfor iPar=1:nParToDo
% for iPar=4:4
    %%---------------------------------------------------------------------
    M=Mpar{iPar};
    M.mod{M.i_mod.ModMembrane}.pm.n_val_min=5;
    M.mod{M.i_mod.ModMembrane}.pm.n_val_max=8;
%     M.mod{M.i_mod.ModMembrane}.var.CanSplitMerge=ones(M.mod{M.i_mod.ModMembrane}.var.n_edg,2);
%     M.mod{M.i_mod.ModMembrane}.pm.remeshExtraCtrl=1000;
%     M.mod{M.i_mod.ModMembrane}.pm.doRemeshExtraCtrl=true;
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
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirFolderPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
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
save([dirRoot filesep 'CLadding' filesep 'progress1.mat']);
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
save([dirRoot filesep 'CLadding' filesep 'progress2.mat']);
disp('done!');