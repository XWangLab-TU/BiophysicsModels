%function [] = ClathrinMediatedEndocytosis(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % ClathrinMediatedEndocytosis performs the dynamics of CALM,
        % Clathrin and membrane during endocytosis
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, membraneFussion, exampleLamellipodia,
        %   exampleEndocytosis, redBloodCell
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------        
% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('dirRoot', @(x) ischar(x));
% ip.addRequired('dirMod', @(x) ischar(x));
% ip.addParameter('nStep', 20000, @isnumeric); %simulation steps for membrane relaxation
% ip.addParameter('nRep', 12, @isnumeric); %repeat number
% ip.addParameter('xyzLim', [-20 20;-20 20;-20 20], @isnumeric); %mesh range
% ip.addParameter('xyzLimPlot', [-10 10;-10 10;-15 15], @isnumeric); %mesh range
% ip.addParameter('viewAng', [0 0], @isnumeric); %figure view angle
% ip.addParameter('idInit', [], @isnumeric); %which initial condition to use
% ip.addParameter('Tcutoff', 0.00000005, @isnumeric); %time scale of plotting
% ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/CME2';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
% nStep=2000;nRep=12;xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-10 10;-10 10;-15 15];idInit=1;
%--------------------------------------------------------------------------
% nStep=ip.Results.nStep;
% nRep=ip.Results.nRep;
% xyzLim=ip.Results.xyzLim;
% viewAng=ip.Results.viewAng;
% xyzLimPlot=ip.Results.xyzLimPlot;
% idInit=ip.Results.idInit;
% Tcutoff=ip.Results.Tcutoff;%0.00000005
%--------------------------------------------------------------------------
%This code is currently under development in script mode 
%%
% adding objects to path
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
%  load('/endosome/work/bioinformatics/s171152/data/CMEmodel/matlab_epsin.mat')
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/matlab3.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_adjust_k3.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_test_dyn2.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/dyn_k1=20_30.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/dyn_k1=10_20.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/dyn_k1=30_40.mat');
% load('/endosome/work/bioinformatics/s171152/data/CMEmodel/init_k1=20_30_40.mat');
%loading finished initial complex by ClathrinMediatedEndocytosisInit
dirRoot='/endosome/work/bioinformatics/s171152/data/CMEmodel';
iClathrin=20; %having 20 triskelia
dirFolder=[dirRoot filesep 'CLadding' filesep num2str(iClathrin)]; 
recAdd=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],'deleteFiles',false);    
nPar=numel(recAdd.dir_all);
dirPar=cell(nPar,1);
Mpar=cell(nPar,1);
    for iPar=1:nPar
        dirPar{iPar}=recAdd.dir_all{iPar};
        dirTem=dir(recAdd.dir_all{iPar}.dir_data);
        S=load([dirTem(end).folder filesep dirTem(end).name]);
        Mpar{iPar}=S.M;
        clear S;
    end
%--------------------------------------------------------------------------
% parameter setting
nStep=200000; 
k1=[10,20,30]; Nk1=numel(k1);
% k1=[20,30,40]; Nk1=numel(k1);
adjustK3=false; k3adj=300; 
adjustK2=false; k2adj=150;
addDyn=false; kDyn=[5 1];
% whether load half-finished work
loadDone=true;

if adjustK2==true
    nameCond=['K2=' num2str(k2adj)];
elseif addDyn==true
    nameCond=['Dyn2_kD=' num2str(kDyn(1)) '_' num2str(kDyn(2))];
elseif adjustK3==true
    nameCond=['K3=' num2str(k3adj)];
else
    nameCond='WT';
end
    
dirTem=[dirRootSt3 '_' nameCond '_k1='];
for i=1:Nk1
    dirTem=[dirTem '_' num2str(k1(i))];
end

if loadDone==true
    rec=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',12,'deleteFiles',false);
else
    rec=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',12,'deleteFiles',true);
end
nPar=numel(rec.dir_all);
i_loop_start=cell(nPar,1);
for iPar=1:nPar
i_loop_start{iPar}=0;
end
FnamePar=cell(nPar,1);
for iPar=1:nPar
    if loadDone==true
    %----------------------------------------------------------------------
    dirInfo=dir(rec.dir_all{iPar}.dir_data);
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=0;
    for iFile=1:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    if nIdx>nIdxMax
        nIdxMax=nIdx;
        iFileRead=iFile;
    end
    end    
    S=load([dirInfo(iFileRead).folder filesep dirInfo(iFileRead).name]);
    i_loop_start{iPar}=nIdxMax;
    Mpar{iPar}=S.M;
    end
    %----------------------------------------------------------------------
    Mpar{iPar}.TypForce.pm.k_ModClathrin_ModFreeParticle=200;
%     Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=20*3; %100
%     Mpar{iPar}.TypForce.pm.k_ModClathrin(2)=200.; %1000
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=50000+10000*iPar;
%     Mpar{iPar}.TypForce.otherInfo.ModClathrin_ModFreeParticle.idC_idFP(1:3:end,:)=[];
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.P=0.;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=0.;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=100+0*iPar; %100
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=200000; %500000
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %0.2
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
    Mpar{iPar}.TypForce.int_stored.ModMembrane=[];
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModFreeParticle}.prop = {'Particle'};
%--------------------------------------------------------------------------
    if addDyn==false
    FnamePar{iPar}={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    else
    FnamePar{iPar}={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModFreeParticle_ModMembrane'};
    end
%--------------------------------------------------------------------------
end
for iPar=1:nPar          
    dirPar{iPar}=rec.dir_all{iPar};
end
for iK1=1:Nk1
    iParStart=(iK1-1)*nPar/Nk1+1;
    iParEnd=(iK1)*nPar/Nk1;
    for iPar=iParStart:iParEnd
    Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=k1(iK1);
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=500000/200*Mpar{iPar}.TypForce.pm.k_ModClathrin(1);
    end
end
for iPar=1:nPar
    if adjustK2==true
        Mpar{iPar}.TypForce.pm.k_ModClathrin(2)=k2adj;
    end
    if adjustK3==true
        Mpar{iPar}.TypForce.pm.k_ModClathrin_ModFreeParticle=k3adj;
    end
    if addDyn==true
        Mpar{iPar}.TypForce.pm.k_ModFreeParticle_ModMembrane=kDyn;
    end
end

fprintf('relaxing. ..\n');
parfor iPar=1:nPar
% for iPar=nPar:nPar
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname=FnamePar{iPar};
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModFreeParticle_ModMembrane'};
%     Fname={'ModMembrane'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=i_loop_start{iPar};
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',0.00000005,'plot_or_not',false,'sMax', 0.5);
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false,'iStage',2); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    loop_done_pre{iPar}=i_loop_done;
    Mpar{iPar}=M;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------