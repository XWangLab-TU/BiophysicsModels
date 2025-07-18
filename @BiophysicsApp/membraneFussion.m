function [dataStored] = membraneFussion(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % membraneFussion performs the two vesicle fusing into one vesicle case study
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, exampleLamellipodia,
        %   exampleFilopodia, exampleEndocytosis, redBloodCell
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', 50000, @isnumeric); %simulation steps for membrane relaxation
ip.addParameter('nRep', 12, @isnumeric); %repeat number
ip.addParameter('xyzLim', [-7 13;-10 10;-10 10], @isnumeric); %figure axis range
ip.addParameter('viewAng', [15 10], @isnumeric); %figure view angle
ip.addParameter('CSremesh', false, @islogical); %to test CS-based remeshing
ip.addParameter('CSremeshAnalysisOnly', false, @islogical);
ip.addParameter('GaussianCurv', false, @islogical);
ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig3_fussion';
% dirMod='/home2/s171152/codes/matlab/mine/git/module/module';addpath(dirMod);
%--------------------------------------------------------------------------
nStep=ip.Results.nStep;
nRep=ip.Results.nRep;
xyzLim=ip.Results.xyzLim;
viewAng=ip.Results.viewAng;
CSremesh=ip.Results.CSremesh;
CSremeshAnalysisOnly=ip.Results.CSremeshAnalysisOnly;
GaussianCurv=ip.Results.GaussianCurv;
%--------------------------------------------------------------------------
if CSremesh==true
    nGroup=2;
    dirRoot1=[dirRoot filesep 'CS'];
    dirRoot2=[dirRoot filesep 'ORG'];
elseif GaussianCurv==true
    nGroup=2;
    dirRoot1=[dirRoot filesep 'GC']; %Gaussian curvature
    dirRoot2=[dirRoot filesep 'ORG'];
else
    nGroup=1;
end
%==========================================================================
if (CSremeshAnalysisOnly==false)
for iGroup=1:nGroup
%--------------------------------------------------------------------------
if CSremesh==false
    dirRootTem=dirRoot;
elseif (CSremesh==true) && (iGroup==1)
    dirRootTem=dirRoot1;
elseif (CSremesh==true) && (iGroup==2) 
    dirRootTem=dirRoot2;
end
rec=ComRecord([dirRootTem filesep 'data'],...
              [dirRootTem filesep 'fig'],...
              [dirRootTem filesep 'init'],'nRep',nRep,'deleteFiles',true);
nPar=numel(rec.dir_all);    
%--------------------------------------------------------------------------
%Compute membrane fussion
%--------------------------------------------------------------------------
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); %baseline length 1000;
m=ModMembrane(true,2,0,'unit',u);
if (CSremesh==true) && (iGroup==1)
    m.pm.remeshScheme=1; %CS-based remesh
end
if (GaussianCurv==true) && (iGroup==1)
    m.pm.GaussianCurv=true; %consider Gaussian curvature
end
m.pm.n_val_max=10; %maximal valence reaches 10 right after the fusion process
m2=m; %copy the vesicle but shift it to the right of the 1st vesicle
m2.var.coord(:,1)=m2.var.coord(:,1)+(max(m.var.coord(:,1))-min(m2.var.coord(:,1)));
%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],m);
%docking the 2 vesicles
[M] = M.mod{M.i_mod.ModMembrane}.docking(m,m2,M);
    %adjust parmeters for this particular application
    %------------------------------------------
    M.mod{M.i_mod.ModMembrane}.pm.n_val_max=10; %maximal valence reaches 10 right after the fusion process
    M.mod{M.i_mod.ModMembrane}.pm.k_V=4; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_A=8; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_a=0; %12
    M.mod{M.i_mod.ModMembrane}.pm.P=0.; %0
    M.mod{M.i_mod.ModMembrane}.pm.dt=0.0025;
    M.mod{M.i_mod.ModMembrane}.pm.k_c=5; %5
    M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.2; %0.2
    M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.001;
    M.mod{M.i_mod.ModMembrane}.pm.kBT=0.0; 
    M.mod{M.i_mod.ModMembrane}.pm.V0=312; %volume of 2 spheres
    M.mod{M.i_mod.ModMembrane}.pm.A0=280; %area of 2 spheres
    M.mod{M.i_mod.ModMembrane}.pm.nAVmean=4;
    %------------------------------------------
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
dataStored=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
    dirPar{iPar}=rec.dir_all{iPar};
    dataStored{iPar}=[];
end
%--------------------------------------------------------------------------
%dynamics
fprintf('computing vesicle fusion morphology...\n');
%plot initial state
parfor iPar=1:nPar
    if CSremesh==false
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,0,'rec_fig_only',false,'update_mod_only',false);
    end
end
parfor iPar=1:nPar
% for iPar=11:11    
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModMembrane'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    %stop simulation at every Tcutoff to save results
    Tcutoff=M.mod{M.i_mod.ModMembrane}.pm.dt*0.2;
    while (i_loop_done<nStep)
        coordSave=[];
        edgSave=[];
        if CSremesh==true
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
        coordSave=M.mod{M.i_mod.ModMembrane}.var.coord;
        edgSave=M.mod{M.i_mod.ModMembrane}.var.edge_all;
        end
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'sMax',10);
        dataStored{iPar}=[dataStored{iPar};CutOff.dataStored];
        if (CSremesh==true)
            if isempty(CutOff.dataStored)
                nEdg=0;
            else
            nEdgTot=size(CutOff.dataStored,1);
            for iEdg=1:nEdgTot
                if CutOff.dataStored(iEdg,1)==CutOff.dataStored(1,1)
                    nEdg=iEdg;
                else
                    break;
                end
            end
            idEdg=CutOff.dataStored(1:nEdg,2);
            end
            figure(fig);hold on;
            for i=1:nEdg
                plot3(coordSave(edgSave(idEdg(i),:),1),coordSave(edgSave(idEdg(i),:),2),coordSave(edgSave(idEdg(i),:),3),...
                      'color',[0 1 0], 'linewidth',2);
                hold on;
            end
            ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
        end
        i_loop_done=i_loop_done+CutOff.nS;
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
        if CSremesh==false
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
        end
        %--------------------------------------------------------------------------
    end
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('finishded at i_par %d \n',iPar);
end
if (CSremesh==true) && (iGroup==1)
    save([dirRoot1 filesep 'DataRes.mat'],'dataStored');
elseif (CSremesh==true) && (iGroup==2)
    save([dirRoot2 filesep 'DataRes.mat'],'dataStored');
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end
end
%%
if CSremesh==true
S1=load([dirRoot1 filesep 'DataRes']);
S2=load([dirRoot2 filesep 'DataRes']);
dataComp=struct('nRemesh',cell(2,1),'coordVar',cell(2,1),'H',cell(2,1));

for i=1:nRep
    dataComp(1).nRemesh=[dataComp(1).nRemesh;size(S1.dataStored{i},1)];
    dataComp(2).nRemesh=[dataComp(2).nRemesh;size(S2.dataStored{i},1)];
    coord1=0.5*(S1.dataStored{i}(:,5:7)+S1.dataStored{i}(:,8:10));
    coord2=0.5*(S2.dataStored{i}(:,5:7)+S2.dataStored{i}(:,8:10));
    dataComp(1).coordVar=[dataComp(1).coordVar;mean(sqrt(sum(coord1.^2,2)))];
    dataComp(2).coordVar=[dataComp(2).coordVar;mean(sqrt(sum(coord2.^2,2)))];
    H1=0.5*(S1.dataStored{i}(:,3)+S1.dataStored{i}(:,4));
    H2=0.5*(S2.dataStored{i}(:,3)+S2.dataStored{i}(:,4));
    dataComp(1).H=[dataComp(1).H;mean(H1)];
    dataComp(2).H=[dataComp(2).H;mean(H2)];
end
%%
colBar=[0 1 1; 1 1 0; 1 0 1];
colTem=[0,0,0]; colTem2=[1,0,0];
   markerStyle={'o','s','p'};
   nPlot=2;
   for iPlot=1:nPlot
        fig=figure;
        if iPlot==1
            dataToPlot=log10(dataComp(1).nRemesh);
        elseif iPlot==2
            dataToPlot=(dataComp(1).coordVar);
        else
            dataToPlot=(dataComp(1).H);
        end
        H=ComPlot.notBoxPlot(dataToPlot,0,'jitter',2); hold on;
        set([H.data],'Marker',markerStyle{1},'MarkerFaceColor',colTem,'markerEdgeColor',colTem,'MarkerSize',3);
        set([H.mu],'color',colTem2);
        colP=colBar(1,:);          
        set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
        set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
        if iPlot==1
            dataToPlot=log10(dataComp(2).nRemesh);
        elseif iPlot==2
            dataToPlot=(dataComp(2).coordVar);
        else
            dataToPlot=(dataComp(2).H);
        end
        H=ComPlot.notBoxPlot(dataToPlot,3,'jitter',2); hold on;
        set([H.data],'Marker',markerStyle{2},'MarkerFaceColor',colTem,'markerEdgeColor',colTem,'MarkerSize',3);
        set([H.mu],'color',colTem2);
        colP=colBar(2,:);          
        set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
        set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
        xticks([0 3]);
        ticklabels={'Geo', 'FE'};
        xticklabels(ticklabels);
        if iPlot==1
            ylabel('$log_{10}N_{remesh}$','interpreter','latex');
        elseif iPlot==2
            ylabel('$\sigma (\vec{|r|})$','interpreter','latex');
        else
            ylabel('$\sigma (\vec{f})$','interpreter','latex');
        end
   end
savefig(fig,[dirRoot filesep 'AnalysisRes.fig']);
end
%==========================================================================
end