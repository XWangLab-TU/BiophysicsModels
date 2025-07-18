function [] = tetherPull(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % tetherPull performs the tether pulling case study
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also membraneFussion, exampleLamellipodia,
        %   exampleFilopodia, exampleEndocytosis, redBloodCell
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('kDiffAlt', [0.025 0.05 0.075 0.1 0.125 0.15], @isnumeric); %various diffusion constant for tension propagation
ip.addParameter('r_threAlt', 0.1, @isnumeric); %a small gap to allow diffusion to pass from left tether to right tether
ip.addParameter('iMechanism', [1 2 3], @isnumeric); %1: instant tension propagation; 2: diffusion; 3: diffusion and extra membrane
ip.addParameter('reloadInit', true, @islogical); %true: don't compute initial tether
ip.addParameter('idInitNormal', [], @isnumeric); %which initial tether to use
ip.addParameter('reloadInitExtraMem', true, @islogical); %true: don't compute tether with extra membrane
ip.addParameter('idInitExtraMem', 10, @isnumeric); %which tether with extra membrane to use, 10 being the most extra; use [] to see all
ip.addParameter('nRepInitNormal', 10, @isnumeric); %how many initial tether to compute
ip.addParameter('nRepInitExtraMem', 10, @isnumeric); %how many tether with extra membrane to compute
ip.addParameter('nRepPull', 8, @isnumeric); %how many tether pulling per each parameter setting to compute
ip.addParameter('stepInit', 40000, @isnumeric); %simulation steps for inital tether relaxation
ip.addParameter('stepPull', 20000, @isnumeric); %simulation steps for pulling tether relaxation
ip.addParameter('reloadComputed', true, @islogical); %avoid repeating previous computations
ip.addParameter('reloadAnalysis', true, @islogical); %avoid repeating previous analysis
ip.addParameter('iStartMeasure', 6, @isnumeric); %defines position to start measuring tether radius
ip.addParameter('iEndMeasureL', 30, @isnumeric); %defines position to end measuring left tether radius
ip.addParameter('iEndMeasureR', 16, @isnumeric); %defines position to end measuring right tether radius
ip.addParameter('plotMeasuredRings', 0, @isnumeric); %plot measured rings: 0-no, 1-plot first one, 2-everyone
ip.addParameter('rMeasureMethod', 0, @isnumeric); %0-mean, 1-median, 2-max of multiple rings around tethers
ip.addParameter('scaleBarPlot', 100, @isnumeric); %adjust bar width vs plot width
ip.addParameter('frameMax', 30, @isnumeric);
ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig5';
% dirMod='/home2/s171152/codes/matlab/mine/git/module/module';addpath(dirMod);
%--------------------------------------------------------------------------
kDiffAlt=ip.Results.kDiffAlt;
r_threAlt=ip.Results.r_threAlt;
iMechanism=ip.Results.iMechanism;
reloadInit=ip.Results.reloadInit;
idInitNormal=ip.Results.idInitNormal;
reloadInitExtraMem=ip.Results.reloadInitExtraMem;
idInitExtraMem=ip.Results.idInitExtraMem;
nRepInitNormal=ip.Results.nRepInitNormal;
nRepInitExtraMem=ip.Results.nRepInitExtraMem;
stepInit=ip.Results.stepInit;
stepPull=ip.Results.stepPull;
reloadComputed=ip.Results.reloadComputed;
reloadAnalysis=ip.Results.reloadAnalysis;
nRepPull=ip.Results.nRepPull;
iStartMeasure=ip.Results.iStartMeasure;
iEndMeasureL=ip.Results.iEndMeasureL;
iEndMeasureR=ip.Results.iEndMeasureR;
plotMeasuredRings=ip.Results.plotMeasuredRings;
rMeasureMethod=ip.Results.rMeasureMethod;
scaleBarPlot=ip.Results.scaleBarPlot;
frameMax=ip.Results.frameMax;
%==========================================================================
% compute inital 2 tethers on +-45 degree of a sphere, also with extra
% membrane
%==========================================================================
%--------------------------------------------------------------------------
%setup path using @ComRecord
if reloadInit==true
recInitNormal=ComRecord([dirRoot filesep 'normalTether' filesep 'data'],...
              [dirRoot filesep 'normalTether' filesep 'fig'],...
              [dirRoot filesep 'normalTether' filesep 'init'],'nRep',nRepInitNormal,'deleteFiles',false);
else
recInitNormal=ComRecord([dirRoot filesep 'normalTether' filesep 'data'],...
              [dirRoot filesep 'normalTether' filesep 'fig'],...
              [dirRoot filesep 'normalTether' filesep 'init'],'nRep',nRepInitNormal,'deleteFiles',true);    
end
if reloadInitExtraMem==true
    recInitExtraMem=ComRecord([dirRoot filesep 'extraMemTether' filesep 'data'],...
              [dirRoot filesep 'extraMemTether' filesep 'fig'],...
              [dirRoot filesep 'extraMemTether' filesep 'init'],'nRep',nRepInitExtraMem,'deleteFiles',false);
else
    recInitExtraMem=ComRecord([dirRoot filesep 'extraMemTether' filesep 'data'],...
              [dirRoot filesep 'extraMemTether' filesep 'fig'],...
              [dirRoot filesep 'extraMemTether' filesep 'init'],'nRep',nRepInitExtraMem,'deleteFiles',true);
end
nPar=numel(recInitNormal.dir_all);
%--------------------------------------------------------------------------
%Compute normal tether
%--------------------------------------------------------------------------
if reloadInit==false
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); %baseline length 1000;
m=ModMembrane(true,3,0,'unit',u);
s=ModSubstrate(0.01,13);
%assemble the model with the membrane and substrate (puller and holder external points) objects
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{m,s}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],m,s);
    %adjust parmeters for this particular application
    %------------------------------------------
    M.mod{M.i_mod.ModMembrane}.pm.k_V=4; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_A=0; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_a=12; %12
    M.mod{M.i_mod.ModMembrane}.pm.P=0.; %0
    M.mod{M.i_mod.ModMembrane}.pm.dt=0.005;
    M.mod{M.i_mod.ModMembrane}.pm.k_c=5; %5
    M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.2; %0.2
    M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.001;
    M.mod{M.i_mod.ModMembrane}.pm.kBT=0.0; %0.01
    M.mod{M.i_mod.ModMembrane}.pm.V0=1223; %1223*0.6,9735
    M.mod{M.i_mod.ModMembrane}.pm.A0=553; %553,2204
    M.mod{M.i_mod.ModMembrane}.pm.nAVmean=4;
    M.TypForce.pm.k_ModMembrane_ModSubstrate=5;
    M.mod{M.i_mod.ModMembrane}.pm.kDiff=inf;
    %------------------------------------------
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
%     S=load('/home2/s171152/codes/matlab/mine/git/module/temp/matlab.mat');
%     Mpar{iPar}=S.mod_save;
    Mpar{iPar}=M;
    dirPar{iPar}=recInitNormal.dir_all{iPar};
end
%--------------------------------------------------------------------------
%dynamics
fprintf('computing normal tether\n');
parfor iPar=1:nPar
% for iPar=1:1    
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(stepInit);
    Fname={'ModMembrane';'ModMembrane_ModSubstrate'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    %stop simulation at every Tcutoff to save results
    Tcutoff=M.mod{M.i_mod.ModMembrane}.pm.dt*0.1;
    while (i_loop_done<stepInit)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'sMax',10);
        i_loop_done=i_loop_done+CutOff.nS;
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', [-25 25;-25 25;-25 25],'viewAng',[0 90]);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',true);
        %--------------------------------------------------------------------------
    end
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('finishded at i_par %d \n',iPar);
end
idInitNormal=[];
end
%--------------------------------------------------------------------------
%choose a case from repeats for next steps
if isempty(idInitNormal)
    for iPar=1:nPar
        myDir=dir(recInitNormal.dir_all{iPar}.dir_fig);   
        figTem=open([recInitNormal.dir_all{iPar}.dir_fig filesep myDir(end-1).name]);
        figTem.Name=['init# ' num2str(iPar)];
    end
    prompt='Please type # of inital tether based on prompt figures: ';
    idInitNormal=input(prompt);
    close all;
end
%--------------------------------------------------------------------------
%Compute tether with extra membrane
%--------------------------------------------------------------------------
if reloadInitExtraMem==false
    %load normal tether
    S=load([recInitNormal.dir_all{idInitNormal}.dir_data filesep 'mod_save.mat']);
    M=S.M;
    clear S;
    Mpar=cell(nPar,1);
    dirPar=cell(nPar,1);
    for iPar=1:nPar
        Mpar{iPar}=M;
        dirPar{iPar}=recInitExtraMem.dir_all{iPar};
        %set extra membrane
        %look for vertices on tethers
%         idTem=((M.mod{M.i_mod.ModMembrane}.var.coord(:,1)>6.5) & (M.mod{M.i_mod.ModMembrane}.var.coord(:,2)>6.5)) |...
%               ((M.mod{M.i_mod.ModMembrane}.var.coord(:,1)<-6.5) & (M.mod{M.i_mod.ModMembrane}.var.coord(:,2)>6.5));
%         ComPlot.PlotMod(M);
%         scatter3(M.mod{M.i_mod.ModMembrane}.var.coord(idTem,1),...
%                  M.mod{M.i_mod.ModMembrane}.var.coord(idTem,2),...
%                  M.mod{M.i_mod.ModMembrane}.var.coord(idTem,3),20,'filled');
        %vertices outside of tethers are given higher density for
        %splitting, as extra membrane
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.kDiff=0.05;
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.rDiff=1;
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.xbDiff=1.5;
%         Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.var.dens(~idTem)=Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.var.dens(~idTem)*(1+0.0375*iPar);
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.A0=533+20*iPar;
    end
    %dynamics
    fprintf('computing tether with extra membrane\n');
    parfor iPar=1:nPar
        % for iPar=1:1
        M=Mpar{iPar};
        dyn=dynamics(stepInit);
        Fname={'ModMembrane';'ModMembrane_ModSubstrate'};
        [dyn] = Preparation(dyn,M,Fname);
        i_loop_done=0;
        Tcutoff=M.mod{M.i_mod.ModMembrane}.pm.dt*0.1;
        while (i_loop_done<stepInit)
            [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff);
            i_loop_done=i_loop_done+CutOff.nS;           
            fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
            %plot and save
            fig=ComPlot.PlotMod(M, 'xyzLim', [-25 25;-25 25;-25 25],'viewAng',[0 90]);
            ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',true);
            %--------------------------------------------------------------------------
        end
        Mpar{iPar}=M;
        %makeMovie(rec,'DelayTime', 0.25);
        fprintf('finishded at i_par %d \n',iPar);
    end
    idInitExtraMem=[];
end
%--------------------------------------------------------------------------
%choose a case from repeats for next steps
if isempty(idInitExtraMem)
    for iPar=1:nPar
        myDir=dir(recInitExtraMem.dir_all{iPar}.dir_fig);   
        figTem=open([recInitExtraMem.dir_all{iPar}.dir_fig filesep myDir(end-1).name]);
        figTem.Name=['init# ' num2str(iPar)];
    end
    prompt='Please type # of inital tether based on prompt figures: ';
    idInitExtraMem=input(prompt);
    close all;
end
%==========================================================================
% compute the pulling study case: pulling the ether on the left
%==========================================================================
%scan for two parameters: kDiffAlt (tension diffusion rate) and 
% iMechanism (1: no & 2: diffusion propagation of tension; 3: diffusion propagation with extra tension)
%setup path accordingly
if reloadComputed==true
rec=ComRecord([dirRoot filesep 'data'],[dirRoot filesep 'fig'],[dirRoot filesep 'init'],...
              'nRep',nRepPull,'deleteFiles',false,...
              'paraAlt',struct('kDiffAlt',kDiffAlt,'iMechanism',iMechanism));
else
rec=ComRecord([dirRoot filesep 'data'],[dirRoot filesep 'fig'],[dirRoot filesep 'init'],...
              'nRep',nRepPull,'deleteFiles',true,...
              'paraAlt',struct('kDiffAlt',kDiffAlt,'iMechanism',iMechanism));    
end
nPar=numel(rec.dir_all);
%--------------------------------------------------------------------------
%first, check which conditions need to be computed based on finished result
compOrNot=cell(nPar,1);
if reloadComputed==true
    for iPar=1:nPar
        compOrNot{iPar}=true;
        myDir=dir(rec.dir_all{iPar}.dir_data);
        if numel(myDir)>2
            nFinished=str2double(myDir(end).name(end-8:end-4));
            if nFinished>stepPull
                compOrNot{iPar}=false;
            end
        end
    end
else
    for iPar=1:nPar
        compOrNot{iPar}=true;
    end
end
if reloadComputed==true
    for iPar=1:nPar
        if compOrNot{iPar}==true
           delete([rec.dir_all{iPar}.dir_data filesep '*.*']);
           delete([rec.dir_all{iPar}.dir_fig filesep '*.*']);
        end
    end
end
%--------------------------------------------------------------------------
%read computed tether and setup parameters for 3 diffrent pulling: nothing,
%diffusion barrier and diffusion barrier+extra membrane
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    if rec.dir_all{iPar}.paraVal(2)<3
    S=load([recInitNormal.dir_all{idInitNormal}.dir_data filesep 'mod_save.mat']);
    Mpar{iPar}=S.M;
    clear S;
    else
    S=load([recInitExtraMem.dir_all{idInitExtraMem}.dir_data filesep 'mod_save.mat']);
    Mpar{iPar}=S.M;
    clear S;
    end
    if rec.dir_all{iPar}.paraVal(2)==1
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.rDiff=1;
    else
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.rDiff=r_threAlt;
    end
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.xbDiff=1.5;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.kDiff=rec.dir_all{iPar}.paraVal(1);
    dirPar{iPar}=rec.dir_all{iPar};
end
%--------------------------------------------------------------------------
%dynamics
%%
parfor iPar=1:nPar
    %%
        % for iPar=1:1
        if compOrNot{iPar}==true
        M=Mpar{iPar};
        dyn=dynamics(stepPull);
        Fname={'ModMembrane';'ModMembrane_ModSubstrate'};
        [dyn] = Preparation(dyn,M,Fname);
        i_loop_done=0;
        Tcutoff=M.mod{M.i_mod.ModMembrane}.pm.dt*0.025; %0.025 defines speed of pulling
        while (i_loop_done<stepPull)
            M_save=M;
            %move a external point after every Tcutoff
            M.mod{M.i_mod.ModSubstrate}.var.coord(1,:)=M.mod{M.i_mod.ModSubstrate}.var.coord(1,:)+[-1.5 1.5 0]*0.5; %pulling left tether
            %dynamics
%             CutOff=[];
%             try
            [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff);
%             catch
%                 warning('error found');
%                 fprintf('computing at i_par %d \n',iPar);
%             end
            %--------------------------------------------------------------------------
            %prevent topological defects, if it happens retrieve previous
            %saved result. the defect is when the tether becomes too
            %narrow multiple triangles form a sigularity vertex connecting
            %two pyramids 
            if CutOff.forcedOff==false
            i_loop_done=i_loop_done+CutOff.nS;           
            fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
            %plot add record
            fig=ComPlot.PlotMod(M, 'xyzLim', [-25 25;-25 25;-25 25],'viewAng',[0 90]);
            ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
            else
                M=M_save;
            end
            %--------------------------------------------------------------------------
        end
        Mpar{iPar}=M;
        %makeMovie(rec,'DelayTime', 0.25);
        fprintf('finishded at i_par %d \n',iPar);
        end
end    
%%
%==========================================================================
% analyze the results
%==========================================================================
pathAnalysis=fullfile(rec.dir_all{1,1}.dir_data,'..','..');
nTotComp=nPar;
nTotPval=nTotComp/nRepPull;
nPara=numel(rec.pNameForSearch);
if reloadAnalysis==true   
    S=load([pathAnalysis filesep 'AnalysisRes.mat']);
    res=S.res;
    clear S;
else
    %----------------------------------------------------------------------
    %prepare analysis data structure
    res=cell(nTotPval,1);
    for i=1:nTotPval
       res{i}=struct('rTetherLmean',[],'rTetherRmean',[],...
                     'rTetherLmedian',[],'rTetherRmedian',[],...
                     'rTetherLmax',[],'rTetherRmax',[],...
                     'corrE',[],'ResToRecIdx',[],'nFrame',[],'rdyForAnalysis',[]);
    end
    %----------------------------------------------------------------------
    %collect all different parameter values and point them to their
    %corresponding components in rec
    iChecked=1;
    paraValchecked=nan(nTotPval,nPara);
    paraValchecked(iChecked,:)=rec.dir_all{1}.paraVal;
    res{iChecked}.ResToRecIdx=1;
    paraValtest=nan(1,nPara);
    for iPar=2:nPar
        paraValtest(1,:)=rec.dir_all{iPar}.paraVal;
        if ~ismember(paraValtest,paraValchecked,'rows')
            iChecked=iChecked+1;
            paraValchecked(iChecked,:)=paraValtest;
            res{iChecked}.ResToRecIdx=iPar;
        else
            res{iChecked}.ResToRecIdx=[res{iChecked}.ResToRecIdx;iPar];
        end
    end
    %----------------------------------------------------------------------
    %check if a given folder contains enough frames for analysis
    for iPval=1:nTotPval
        res{iPval}.nFrame=zeros(nRepPull,1);
        res{iPval}.rdyForAnalysis=false(nRepPull,1);
        for iRep=1:nRepPull
        iRec=res{iPval}.ResToRecIdx(iRep);
        myDir=dir(rec.dir_all{iRec}.dir_data);
        res{iPval}.nFrame(iRep)=numel(myDir)-2;
        if res{iPval}.nFrame(iRep)>0
            nFinished=str2double(myDir(end).name(end-8:end-4));
            if nFinished>stepPull
                res{iPval}.rdyForAnalysis(iRep)=true;
            end
        end
        end
    end
    %----------------------------------------------------------------------
    %scan through all legit folders and do the analysis 
    for iPval=1:nTotPval
%     for iPval=1:1
        fprintf('analyzing %d ==> %d\n',iPval,nTotPval)
        fig=[];
        res{iPval}.rTetherLmean=cell(nRepPull,1);
        res{iPval}.rTetherRmean=cell(nRepPull,1);
        res{iPval}.rTetherLmedian=cell(nRepPull,1);
        res{iPval}.rTetherRmedian=cell(nRepPull,1);
        res{iPval}.rTetherLmax=cell(nRepPull,1);
        res{iPval}.rTetherRmax=cell(nRepPull,1);
        for iRep=1:nRepPull
        iRec=res{iPval}.ResToRecIdx(iRep);
        res{iPval}.rTetherLmean{iRep}=nan(res{iPval}.nFrame(iRep),1);
        res{iPval}.rTetherRmean{iRep}=nan(res{iPval}.nFrame(iRep),1);
        res{iPval}.rTetherLmedian{iRep}=nan(res{iPval}.nFrame(iRep),1);
        res{iPval}.rTetherRmedian{iRep}=nan(res{iPval}.nFrame(iRep),1);
        res{iPval}.rTetherLmax{iRep}=nan(res{iPval}.nFrame(iRep),1);
        res{iPval}.rTetherRmax{iRep}=nan(res{iPval}.nFrame(iRep),1);
        if res{iPval}.rdyForAnalysis(iRep)==true
        %------------------------------------------------------------------
        for iFrame=1:res{iPval}.nFrame(iRep)
            %load computed model M
            myDir=dir(rec.dir_all{iRec}.dir_data);
            S=load([myDir(iFrame+2).folder filesep myDir(iFrame+2).name]);
            M=S.M;
            %determine whether to plot
            plotOrnot=false;
            if (plotMeasuredRings==2) || ((plotMeasuredRings==1)&&((iPval==1) && (iRep==1) && (iFrame==res{iPval}.nFrame(iRep))))
                fig=ComPlot.PlotMod(M, 'xyzLim', [-35 15;-10 40;-25 25],'viewAng',[0 90]);
                plotOrnot=true;
            end
            %measure right tether first
            %define 45degree tether direction unit vector n
            LorR=1; %left or right tether
            d=0.5;
            n=[sqrt(2)*d*LorR sqrt(2)*d 0];
            %cover the right tether with various V0 (positions)
            rMeasured=[];
            parfor i=iStartMeasure:iEndMeasureR
                V0=[5*LorR+(i-1)*d*LorR 5+(i-1)*d 0];
                Mtem=M;
                if plotOrnot==true
                [pCircle]=measureRtube(Mtem.mod{Mtem.i_mod.ModMembrane},n,V0,'plot_or_not', true,'fig',fig);
                else
                [pCircle]=measureRtube(Mtem.mod{Mtem.i_mod.ModMembrane},n,V0);    
                end
                rMeasured=[rMeasured;pCircle];
            end
            res{iPval}.rTetherRmean{iRep}(iFrame)=mean(rMeasured(:,3));
            res{iPval}.rTetherRmedian{iRep}(iFrame)=median(rMeasured(:,3));
            res{iPval}.rTetherRmax{iRep}(iFrame)=max(rMeasured(:,3));
            %measure left tether first
            LorR=-1; %left or right tether
            d=0.5;
            n=[sqrt(2)*d*LorR sqrt(2)*d 0];
            rMeasured=[];
            parfor i=iStartMeasure:iEndMeasureL
                V0=[5*LorR+(i-1)*d*LorR 5+(i-1)*d 0];
                Mtem=M;
                if plotOrnot==true
                [pCircle]=measureRtube(Mtem.mod{Mtem.i_mod.ModMembrane},n,V0,'plot_or_not', true,'fig',fig);
                else
                [pCircle]=measureRtube(Mtem.mod{Mtem.i_mod.ModMembrane},n,V0);    
                end
                rMeasured=[rMeasured;pCircle];
            end
            res{iPval}.rTetherLmean{iRep}(iFrame)=mean(rMeasured(:,3));
            res{iPval}.rTetherLmedian{iRep}(iFrame)=median(rMeasured(:,3));
            res{iPval}.rTetherLmax{iRep}(iFrame)=max(rMeasured(:,3));
        end
        %------------------------------------------------------------------
        end
        end
    end
    clear S;
    %----------------------------------------------------------------------
    save([pathAnalysis filesep 'AnalysisRes.mat'],'res');
end
   dataToPlot=zeros(nTotPval,nRepPull);
   %%
   figure;
%    for iPval=4:4:12
   col=[0 1 1; 1 1 0; 1 0 1];
   for iPval=6:6:18
       for iRep=1:1
       if frameMax < res{iPval}.nFrame(iRep)
           nFrame=frameMax;
       else
           nFrame=res{iPval}.nFrame(iRep);
       end
       x=res{iPval}.rTetherLmean{iRep}(1:nFrame)/res{iPval}.rTetherLmean{iRep}(1);
       y=res{iPval}.rTetherRmean{iRep}(1:nFrame)/res{iPval}.rTetherRmean{iRep}(1);
       p = polyfit(x,y,4);
       yfit =  p(1) * x.^4 + p(2)* x.^3+p(3)* x.^2+p(4)* x+p(5);
       if iPval<=6
           iCol=1;
       elseif iPval<=12
           iCol=2;
       else
           iCol=3;
       end
       hold on;plot(x,y,'.','color',col(iCol,:));hold on;plot(x,yfit,'-','linewidth',2,'color',col(iCol,:));
       hold on;plot(x(end),yfit(end),'o','linewidth',2,'Markersize',8,'Markeredgecolor',col(iCol,:),'MarkerFaceColor','black')
       end
   end
   xlim([0.6 1.1]);ylim([0.6 1.1]);
   savefig(gcf,[pathAnalysis filesep 'example.fig']);
   %%
   for iPval=1:nTotPval
       for iRep=1:nRepPull
           if frameMax < res{iPval}.nFrame(iRep)
               nFrame=frameMax;
           else
               nFrame=res{iPval}.nFrame(iRep);
           end
           if rMeasureMethod==0
                x=res{iPval}.rTetherLmean{iRep}(1:nFrame)/res{iPval}.rTetherLmean{iRep}(1);
                y=res{iPval}.rTetherRmean{iRep}(1:nFrame)/res{iPval}.rTetherRmean{iRep}(1);
           elseif rMeasureMethod==1
                x=res{iPval}.rTetherLmedian{iRep}(1:nFrame)/res{iPval}.rTetherLmedian{iRep}(1);
                y=res{iPval}.rTetherRmedian{iRep}(1:nFrame)/res{iPval}.rTetherRmedian{iRep}(1);
           elseif rMeasureMethod==2
                x=res{iPval}.rTetherLmax{iRep}(1:nFrame)/res{iPval}.rTetherLmax{iRep}(1);
                y=res{iPval}.rTetherRmax{iRep}(1:nFrame)/res{iPval}.rTetherRmax{iRep}(1);
           else
                warning('wrong rMeasureMethod but mean method used');
                x=res{iPval}.rTetherLmean{iRep}(1:nFrame)/res{iPval}.rTetherLmean{iRep}(1);
                y=res{iPval}.rTetherRmean{iRep}(1:nFrame)/res{iPval}.rTetherRmean{iRep}(1);
           end
%            p = polyfit(x,y,1);
%            yfit =  p(1) * x + p(2);
%            yresid  = y - yfit;
%            SSresid = sum(yresid.^2);
%            SStotal = (length(y)-1) * var(y);
%            rsq = 1 - SSresid/SStotal;
%            dataToPlot(iPval,iRep)=p(1);
             p = polyfit(x,y,4);
             yfit =  p(1) * x.^4 + p(2)* x.^3+p(3)* x.^2+p(4)* x+p(5);
             dataToPlot(iPval,iRep)=yfit(end);
       end
   end
   colBar=[0 1 1; 1 1 0; 1 0 1];
   markerStyle={'o','s','p'};
   fig=figure;
   for iPval=1:nTotPval
        iRec=res{iPval}.ResToRecIdx(1);
        %center bar based on the diffusion constant kDiff (1st parameter scan)
        ctr=rec.dir_all{iRec}.paraVal(1)*scaleBarPlot;
        %different colors for different mechanisms
        if rec.dir_all{iRec}.paraVal(2)==1
            colP=colBar(1,:);
            idxBar=1;
        elseif rec.dir_all{iRec}.paraVal(2)==2
            colP=colBar(2,:);
            idxBar=2;
        else
            colP=colBar(3,:); 
            idxBar=3;
        end
        H=ComPlot.notBoxPlot(dataToPlot(iPval,:),ctr+abs(idxBar-2),'jitter',2); hold on; 
        if rec.dir_all{iRec}.paraVal(2)==1
            colTem=[0,0,0]; colTem2=[1,0,0];
            set([H.data],'Marker',markerStyle{rec.dir_all{iRec}.paraVal(2)},'MarkerFaceColor',colTem,'markerEdgeColor',colTem,'MarkerSize',3);
            set([H.mu],'color',colTem2);
        elseif rec.dir_all{iRec}.paraVal(2)==2
            colTem=[0,0,0]; colTem2=[1,0,0];
            set([H.data],'Marker',markerStyle{rec.dir_all{iRec}.paraVal(2)},'MarkerFaceColor',colTem,'markerEdgeColor',colTem,'MarkerSize',3);
            set([H.mu],'color',colTem2);
        else
            colTem=[0,0,0]; colTem2=[1,0,0];
            set([H.data],'Marker',markerStyle{rec.dir_all{iRec}.paraVal(2)},'MarkerFaceColor',colTem,'markerEdgeColor',colTem,'MarkerSize',3);
            set([H.mu],'color',colTem2);
        end            
            set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
            set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');

   end
   xticks(kDiffAlt*scaleBarPlot);
   nTick=numel(kDiffAlt);
   ticklabels={};
   for i=1:nTick
       ticklabels=[ticklabels;{num2str(kDiffAlt(i))}];
   end
   xticklabels(ticklabels);
   savefig(fig,[pathAnalysis filesep 'AnalysisRes.fig']);
end