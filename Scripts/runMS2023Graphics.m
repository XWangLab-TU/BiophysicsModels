%==========================================================================
% This script adjust the graphic output for the 2023 manuscript. 
% Thie script is supposed to run after finishing runMS2022Examples
% Applications are separated into sections by double %% signs 
% and can be run individually with 'ctrl'+'enter'
% see the documentation of individual application for more detail
%   See also runMS2022Examples
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/25
%==========================================================================
%% add library
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels'; 
dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
addpath(genpath(dirMod));
%% ==========================================================================
% Graphics
%==========================================================================
% plot 2 clathrin triskilia
%% 1.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[0. 0 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0 0]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=0.25;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',facealpha);
xlim([-3 5]); ylim([-4 4]); zlim([-6 2]);
axis off
view([5 0])
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 1.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[0. 1.2 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0.3 -1.3]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=1;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',1);
xlim([-3.5 3.5]); ylim([-3.5 3.5]); zlim([-5.5 1.5]);
axis off
view([30 15])
%% 2.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[-0.2 1 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0 -1]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=0.25;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',facealpha);
xlim([-3 5]); ylim([-4 4]); zlim([-6 2]);
axis off
view([0 50])
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 3.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
m=ModMembrane(true,4,0,'unit',u);
xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-13.5 13.5;-13.5 13.5;-13.5 13.5];
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{},{}};{{}}});
M = model(int_info,dirMod,u,xyzLim,m);
facealpha=1;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
axis off
% xlabel('x (10 nm)');ylabel('y (10 nm)');zlabel('z (10 nm)');
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 3.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
m=ModMembrane(true,4,0,'unit',u);
c2=ModClathrin('unit',u);
c2.var.a(ic,:)=c2.var.a(ic,:)+[-0.2 0.5 0.2];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[0 2 18]; 
[c2.var] = setVar(c2,true);  
xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-15 15;-15 15;-15 15];
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{},{}};{{}}});
M = model(int_info,dirMod,u,xyzLim,m);
facealpha=1;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
plot(c2,'f',fig,'col',[0 0 1],'facealpha',facealpha);
axis off
% xlabel('x (10 nm)');ylabel('y (10 nm)');zlabel('z (10 nm)');
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% plot initial membrane morphologies: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep04/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep11/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep05/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
%% plot initial membrane morphologies showing control points: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/initControlPoints.mat');
Mpar=S.Mpar;
for iPar=1:4:12
    M=Mpar{iPar};
    ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModSubstrate]);
    xlim([-8 8]);ylim([-4 12]);zlim([-14 -6]); axis off; view([0 10]);
end
%% plot initial clathrin-membrane complex: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_10_20/data/rep01/mod_Stg_002Fr_00923');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep01/mod_Stg_002Fr_00018');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_10_20/data/rep01/mod_Stg_002Fr_00754');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
%% plot complex relax: WT
xyzLimPlot=[-8 8;-8 8;-18 -2];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_00927');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_78129');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_200534');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot complex relax: CALM-
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_00937');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_78396');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_196132');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot complex relax: Epsin1-
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_00719');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_75777');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_191895');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

%% phase diagram
% dirAddon='simp';
dirAddon='simp_dyn';
dirMid='020'; nFolder=5; nPara=2; nRep=6;
nSet=nFolder*nPara*nRep;
dirSet=cell(nSet,1);
dataToPlot=struct('K',zeros(nSet,1),'T',zeros(nSet,1));

for iFolder=1:nFolder
    for iPara=1:nPara
        for iRep=1:nRep
            iTem=(iFolder-1)*nPara*nRep+(iPara-1)*nRep+iRep;
            dirSet{iTem}=[dirRoot filesep 'M' num2str(iFolder) dirAddon filesep dirMid filesep 'data' ...
                          filesep 'p' num2str(iPara) '_rep' num2str(iRep)];
        end
    end
end

for iSet=1:nSet
    dirTem=dir(dirSet{iSet});
    S=load([dirTem(end).folder filesep dirTem(end).name]);
    dataToPlot(iSet).K=S.M.TypForce.pm.k_ModClathrin(1);
    dataToPlot(iSet).T= getTightness(S.M.mod{S.M.i_mod.ModClathrin},S.M);
end
clear S;

KT=zeros(nSet,1);
for iSet=1:nSet
    KT(iSet,1)=dataToPlot(iSet).K;
    KT(iSet,2)=dataToPlot(iSet).T;
end

fig=figure;
nBar=nFolder*nPara;
Xscale=4;
colP1=[0 0 1];
colP2=[0 1 0];
BarWadd=0.3;
Njitter=3;
for iFolder=1:nFolder
    for iPara=1:nPara
        iTem=(iFolder-1)*nPara*nRep+(iPara-1)*nRep+1:(iFolder-1)*nPara*nRep+(iPara-1)*nRep+nRep;   
        x=KT(iTem,1);
        y=KT(iTem,2);
        y1=y(y>34);
        y2=y(y<34);
        if ~isempty(y1)
            H=ComPlot.notBoxPlot(y1,x(1),'jitter',Njitter); hold on;
            VertDef=H.sdPtch.Vertices;
            VertNew=VertDef;
            dX=VertDef(2,1)-VertDef(1,1);
            VertNew(1,1)=VertDef(1,1)-BarWadd*dX*0.5;
            VertNew(2,1)=VertDef(2,1)+BarWadd*dX*0.5;
            VertNew(3,1)=VertDef(3,1)+BarWadd*dX*0.5;
            VertNew(4,1)=VertDef(4,1)-BarWadd*dX*0.5;
            set(H.sdPtch,'FaceColor',colP1,'EdgeColor','none','Vertices',VertNew);
            set(H.semPtch,'FaceColor','none','EdgeColor','none');
            set(H.mu,'XData',[VertNew(1,1) VertNew(2,1)],'LineWidth',3);
        end
        if ~isempty(y2)
            H=ComPlot.notBoxPlot(y2,x(1),'jitter',Njitter); hold on;
            VertDef=H.sdPtch.Vertices;
            VertNew=VertDef;
            dX=VertDef(2,1)-VertDef(1,1);
            VertNew(1,1)=VertDef(1,1)-BarWadd*dX*0.5;
            VertNew(2,1)=VertDef(2,1)+BarWadd*dX*0.5;
            VertNew(3,1)=VertDef(3,1)+BarWadd*dX*0.5;
            VertNew(4,1)=VertDef(4,1)-BarWadd*dX*0.5;
            set(H.sdPtch,'FaceColor',colP2,'EdgeColor','none','Vertices',VertNew);
            set(H.semPtch,'FaceColor','none','EdgeColor','none');
            set(H.mu,'XData',[VertNew(1,1) VertNew(2,1)],'LineWidth',3);
        end
    end
end
ylim([15 45]);
% xticks((10:10:100)/Xscale);
%    nTick=numel(xticks);
%    ticklabels={};
%    for i=1:nTick
%        ticklabels=[ticklabels;{num2str(i*10)}];
%    end
% xticklabels(ticklabels);
% xlabel('k_1');
% ylabel('TF');
% xlim([0 110/Xscale])
%% membrane sensitivity
dirAddon='simp_pdf';
% dirAddon='simp_dyn_weak';
% dirAddon='simp_K3_hf';
% dirAddon='simp_K3_db';
dirMid='020'; nFolder=3; nPara=2; nRep=6;
nSet=nFolder*nPara*nRep;
dirSet=cell(nSet,1);
dataToPlot=cell(nSet,1);
steady_state=false;

for iSet=1:nSet
    dataToPlot{iSet}=struct('K',[],'Mmean',[],'Mmedian',[],'Mmax',[],'C',[],'Cleg',[],'T',[]);
end

for iFolder=1:nFolder
    for iPara=1:nPara
        for iRep=1:nRep
            iTem=(iFolder-1)*nPara*nRep+(iPara-1)*nRep+iRep;
            dirSet{iTem}=[dirRoot filesep 'M' num2str(iFolder) dirAddon filesep dirMid filesep 'data' ...
                          filesep 'p' num2str(iPara) '_rep' num2str(iRep)];
        end
    end
end

parfor iSet=1:nSet
    disp(iSet);
    disp(nSet);
    dirTem=dir(dirSet{iSet});
    nFile=size(dirTem,1);
    if steady_state==true
        iFileStart=nFile;
    else
        iFileStart=3;
    end
    %%
    dataToPlot{iSet}.K=zeros(nFile,1);
    dataToPlot{iSet}.T=zeros(nFile,1);
    dataToPlot{iSet}.Mmean=zeros(nFile,1);
    dataToPlot{iSet}.Mmedian=zeros(nFile,1);
    dataToPlot{iSet}.Mmax=zeros(nFile,1);
    dataToPlot{iSet}.C=zeros(nFile,1);
    dataToPlot{iSet}.Cleg=zeros(nFile,1);
    for iFile=iFileStart:nFile
%         disp(iFile);
    S=load([dirTem(iFile).folder filesep dirTem(iFile).name]);
    dataToPlot{iSet}.K(iFile)=S.M.TypForce.pm.k_ModClathrin(1);
    dataToPlot{iSet}.T(iFile)= getTightness(S.M.mod{S.M.i_mod.ModClathrin},S.M);
     M2=getMemWithCoat(S.M);
%     f=M2.TypForce;
%     [f,m,Vtot] = ModMembrane(f,M2);
    f=S.M.TypForce;
    [f,m,Vtot] = ModMembrane(f,S.M);
    [f,V_tot,~] = ModClathrin(f,S.M);
    nLeg=0;
    for i=1:S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord
        nLeg=nLeg+sum(sum(S.M.mod{S.M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
    end
    nLeg=nLeg*0.5;
%     dataToPlot(iSet).Fmem(iFile)= norm(f.int_comp.ModMembrane)/m.var.n_coord;
%     dataToPlot(iSet).Fmem(iFile)= Vtot/m.var.n_coord;
     idx=unique(M2.mod{M2.i_mod.ModMembrane}.var.face_unq);
     dataToPlot{iSet}.Mmean(iFile)=mean(m.var.f.kH(idx).*m.var.f.kH(idx));
     dataToPlot{iSet}.Mmedian(iFile)=median(m.var.f.kH(idx).*m.var.f.kH(idx));
     dataToPlot{iSet}.Mmax(iFile)=max(m.var.f.kH(idx).*m.var.f.kH(idx));
     dataToPlot{iSet}.C(iFile)=V_tot;
     dataToPlot{iSet}.Cleg(iFile)=nLeg;
    end
end
clear S;
%%
nClose=0;
nOpen=0;
for iSet=1:1:nSet
iStart=3;
if dataToPlot{iSet}.T(end)<33
    col=[1 0 0]; nOpen=nOpen+numel(dataToPlot{iSet}.C);
else
    col=[0 0 1]; nClose=nClose+numel(dataToPlot{iSet}.C);
end
scatter((dataToPlot{iSet}.C(3)-dataToPlot{iSet}.C(iStart:end))./dataToPlot{iSet}.Cleg(iStart:end),...
     dataToPlot{iSet}.Mmedian(iStart:end),'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.2,'MarkerFaceColor',col);
hold on;
end
%% histogram
histData=[];
for iSet=1:nSet
    histData=[histData;dataToPlot{iSet}.T(end)];
end
%%
histogram(histData, 'BinWidth', 1,...
'DisplayStyle','stairs');
%% clathrin gliding
xyzLimPlot=[-15 15; -15 15; -28 2];
viewAng=[0 90];
rMaxFromCtr=10;
rMinFromCtr=5;
Tcutoff=5.0000e-08;
D=dir([dirRoot filesep 'M1/008/data/p10_rep1']);
S=load([D(end).folder filesep D(end).name]);
M=S.M;

muSave=M.mod{M.i_mod.ModMembrane}.pm.mu;
M.mod{M.i_mod.ModMembrane}.pm.mu=0.01; %temporarily immobilize the membrane
%----------------------------------------------------------------------
[M,~] = add(M.mod{M.i_mod.ModClathrin},M,[0 0 -3.5],'rMaxFromCtr',rMaxFromCtr,'rMinFromCtr',rMinFromCtr);
fig1=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
           [M.i_mod.ModClathrin,M.i_mod.ModFreeParticle,M.i_mod.ModMembrane]);
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
%     fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
end
%----------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.pm.mu=muSave;
fig2=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
           [M.i_mod.ModClathrin,M.i_mod.ModFreeParticle,M.i_mod.ModMembrane]);
clear S; 
%% ext potential
xyzLimPlot=[-17 17; -17 17; -17 17];
viewAng=[0 0];
D=dir([dirRoot filesep 'M1/002/data/p01_rep1']);
S=load([D(end).folder filesep D(end).name]);
M=S.M;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
           [M.i_mod.ModMembrane]); hold on;
b=-3.25;c=4.25;d=1;a=5;
z=(-10:0.01:10);     
y=-ones(size(z))*16;
plot3(exp(-a*(z-b).^d)+exp(a*(z-c).^d),y,z,'color','g','linewidth',2);
xlim(xyzLimPlot(1,:)); xlim(xyzLimPlot(2,:)); xlim(xyzLimPlot(3,:));
axis off;
clear S;
%% plotting clathrin one by one
dirRoot='/endosome/work/bioinformatics/s171152/data/CMEmodel/CLadding_1';
Para=sprintf('%02d',9); %19 and 9
dirGIF=[dirRoot filesep 'gif'];
mkdir(dirGIF);
zLim=[-15 15];
iCL=1;
CL=sprintf('%03d',iCL);
dirTem=[dirRoot filesep CL filesep 'fig' filesep 'p' Para '_rep1'];
dirData=dir(dirTem);
fig=openfig([dirData(end-1).folder filesep dirData(end-1).name]);
zlim(zLim);
savefig(fig,[dirGIF filesep CL '.fig']);
close(fig);
for iCL=2:23
    CL=sprintf('%03d',iCL);
    dirTem=[dirRoot filesep CL filesep 'fig' filesep 'p' Para '_rep1'];
    dirData=dir(dirTem);
    fig=openfig([dirData(end-1).folder filesep dirData(end-1).name]);
    zlim(zLim);
    savefig(fig,[dirGIF filesep CL '_1' '.fig']);
    close(fig);
    fig=openfig([dirData(3).folder filesep dirData(3).name]);
    zlim(zLim);
    savefig(fig,[dirGIF filesep CL '_0' '.fig']);
    close(fig);
end
ComRecord.makeStackMov(dirGIF,'DelayTime',0.2,'xyzLim',[[-15 15];[-5 25];[-15 15]],'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%--------------------------------------------------------------------------
function M2=getMemWithCoat(M)
       m=M.mod{M.i_mod.ModMembrane};
       AP2=M.mod{M.i_mod.ModFreeParticle};
       nAP2=size(AP2.var.coord,1);
       iM_AP2=[];
       xyzLimPlot(3,2)=max(M.mod{M.i_mod.ModClathrin}.var.coord(:,3))+0.1;
       xyzLimPlot(3,1)=min(M.mod{M.i_mod.ModMembrane}.var.coord(:,3))-0.1;
       for iAP2=1:nAP2
           D=vecnorm(m.var.coord-AP2.var.coord(iAP2,:),2,2);
           [~,iMin]=min(D);
           iM_AP2=[iM_AP2;iMin];
       end  
       iM_AP2=unique(iM_AP2);
       nAP2=numel(iM_AP2);
       iEdgAll=(1:m.var.n_edg)';
       iEdg=[];
       for iAP2=1:nAP2
           iM=iM_AP2(iAP2);
           id=sum(iM==m.var.edge_all,2);
           iEdg=[iEdg;iEdgAll(id>0)];
       end
       iEdg=unique(iEdg);
       nEdg=numel(iEdg);
       iEdgRing=iEdg;
       for iE=1:nEdg
       [~,~,~,id_ring_edg2] = remeshRing(m,iEdg(iE),'ring_ord', 2,'plot_or_not',false);
          iEdgRing=[iEdgRing;id_ring_edg2];
       end
       iEdgRing=unique(iEdgRing);
       iM=m.var.edge_all(iEdgRing,:);
       iM=unique(iM);
       nM=numel(iM);
       nFace=size(m.var.face_unq,1);
       iFaceAll=(1:nFace)';
       iFace=[];
       for i=1:nM
          id=abs(iM(i)-m.var.face_unq);
          id=min(id,[],2);
          id(id>0)=-1;
          id=id+1;
          id=logical(id);
          iFace=[iFace;iFaceAll(id,:)];
       end
       m2=m;
       m2.var.face_unq=m.var.face_unq(iFace,:);
       M2=M;
       M2.mod{M2.i_mod.ModMembrane}=m2;
end