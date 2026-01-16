function [M,fig] = clathrinVesicle(dirLoc,varargin)
%--------------------------------------------------------------------------
        % Author: Xinxin Wang
        % email: wangxinxin8627@gmail.com
        % date: 2026/01/14
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirLoc', @(x) ischar(x));
ip.addParameter('nStep', 200000, @isnumeric);
ip.addParameter('M', [], @isobject);
ip.addParameter('icTot', 20, @isnumeric);
ip.addParameter('k_ModClathrin_ModFreeParticle', 20, @isnumeric);
ip.addParameter('k_ModClathrin', [10 100], @isnumeric);
ip.addParameter('n_init_adp', 2, @isnumeric);
ip.addParameter('meshD', 1, @isnumeric);
ip.addParameter('xyzLimPlot', [-10 10;-10 10;-10 10], @isnumeric);
ip.addParameter('viewAng', [0 0], @isnumeric);
ip.addParameter('r_min_for_link', 2.5, @isnumeric);
ip.parse(dirLoc,varargin{:});
%--------------------------------------------------------------------------
nStep=ip.Results.nStep;
icTot=ip.Results.icTot;
xyzLimPlot=ip.Results.xyzLimPlot;viewAng=ip.Results.viewAng;
if ~isempty(ip.Results.M)
    M=ip.Results.M;
    Vname=cell(2,1);
    Vpm=cell(2,1);
    Vname{1}='VedgMembrane';Vpm{1}=M.mod{M.i_mod.ModMembrane}.pm.Vedg;Vname{2}='VvolMembrane';
    dynLC=dynLattice(1,Vname,Vpm);
    u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
    [lc,~] = ComLattice(M,u);
else
%% add library
% dirLoc='D:\OneDrive\codes\git';
% dirLoc='C:\Users\xiw1834\OneDrive\codes\git';
% dirLoc='/OneDrive/codes/git';
dirMod=[dirLoc filesep 'BiophysicsModels'];
addpath((dirMod));
addpath(([dirLoc filesep 'Internal']));
%cd(dirMod);
%%
% ================================================================ old model
% Mold=load('D:\Matlab_codes\data\mesh\mod_Stg_001Fr_1248628.mat');
% c=Mold.M.mod{Mold.M.i_mod.ModClathrin};
% m=Mold.M.mod{Mold.M.i_mod.ModMembrane};
% p=Mold.M.mod{Mold.M.i_mod.ModFreeParticle};
% f=Mold.M.TypForce;
% ========================================================================
% Computation preparation
%==========================================================================
% unit: 10nm, 1kBT as natural units 
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
% membrane
l0=1;
nVer=1;
m=ModMembrane(true,nVer,0,'unit',u,'l0',l0);
c=ModClathrin('unit',u);
% c.pm.n_min_adp=2;c.pm.n_max_adp=2;
% c.var.coord(:,3)=c.var.coord(:,3)+5;
% later as AP2-adaptor proteins to connect clathrin and membrane
p=ModFreeParticle(0,'unit',u);p.prop = {'Particle','needFollow'};
%introducing various interactions to be included, for example: {m,s} means ctrl points interact with membrane
% int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{c,p},{p,m}};{{}}});
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{c,p}};{{}}});

% D=sqrt(sum((m.var.coord(m.var.edge_all(:,1),:)-m.var.coord(m.var.edge_all(:,2),:)).^2,2)); %histogram(D); m=M.mod{1};
% dr=0.001;
% r2=(m.pm.Vedg.r_1+dr:dr:m.pm.Vedg.r_2-dr)';
% r2=[r2,zeros(size(r2)),zeros(size(r2))];
% r1=zeros(size(r2));
% [V] = potential.VedgMembrane(r1,r2,m.pm.Vedg);
% r=sqrt(sum((r2-r1).^2,2));
% plot(r,V,'.');
% ylim([-0.5 0.1]);
%--------------------------------------------------------------------------
% lattice
% mex timeEvalMex.c; clc;
xyzLim=[-200 200;-200 200;-200 200];
m.pm.nt=1; %4000000
m.pm.k_c=5;

M = model(int_info,dirMod,u,xyzLim,m,c,p);
[lc,~] = ComLattice(M,u);
M = model(int_info,dirMod,u,xyzLim,m,c,p);
Vname=cell(2,1);
Vpm=cell(2,1);
Vname{1}='VedgMembrane';Vpm{1}=M.mod{M.i_mod.ModMembrane}.pm.Vedg;Vname{2}='VvolMembrane';
dynLC=dynLattice(1,Vname,Vpm);

M.TypForce.pm.k_ModClathrin_ModFreeParticle=ip.Results.k_ModClathrin_ModFreeParticle;
M.TypForce.pm.k_ModClathrin=ip.Results.k_ModClathrin;
M.mod{M.i_mod.ModClathrin}.pm.n_init_adp=ip.Results.n_init_adp;
M.mod{M.i_mod.ModClathrin}.pm.r_min_for_link=ip.Results.r_min_for_link;
M.Mesh.d=ip.Results.meshD; 
% inital relax, getting AP2 positions
[lc,M,~] = timeEvalMemCL(dynLC,lc,M);
[M,changed] = docking(M.mod{M.i_mod.ModClathrin},M,[0 0 5],'shiftRange',1);
if ~changed
    error('clathrin init failed!!!');
end
end
%% ========================================================================
% Computation loop of nStep
%==========================================================================
M.mod{M.i_mod.ModMembrane}.pm.nt=nStep;%Msave=M; M=Msave;
for ic=1:icTot
    disp(ic);
%%
% recruit
[M,changed] = recruit(M.mod{M.i_mod.ModClathrin},M,'sortFromCtr', true); 
if ~changed
    error('clathrin recruitment failed!!!');
end
[lc,M,~] = timeEvalMemCL(dynLC,lc,M);
% save('/MATLAB Drive/M.mat',"M");
[M.mod{M.i_mod.ModClathrin},changed] = link(M.mod{M.i_mod.ModClathrin});
if changed
    [lc,M,~] = timeEvalMemCL(dynLC,lc,M);
end
iCfree=sum(sum(sum(M.mod{M.i_mod.ModClathrin}.var.connect==0)));
%%
if (iCfree==0)
    disp('all clathrin legs connected!!!');
    break;
end
end
%% plot
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot', [M.i_mod.ModClathrin,M.i_mod.ModMembrane]);
hold on;
for i=1:M.mod{M.i_mod.ModClathrin}.var.n_coord
    id=M.mod{M.i_mod.ModClathrin}.var.idBond(i,:);id=id(id~=0);
    r=M.mod{M.i_mod.ModFreeParticle}.var.coord(id,:);
    scatter3(r(:,1),r(:,2),r(:,3),40,'filled'); hold on;
end
%==========================================================================
end