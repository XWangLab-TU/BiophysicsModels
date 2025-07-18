function [] = exampleMemCytoskeleton(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % exampleMemCytoskeleton performs the example of Membrane interacting with cytoskeleton
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, membraneFussion, exampleLamellipodia
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2023/03/31
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', 20000, @isnumeric); %simulation steps for membrane relaxation
ip.addParameter('nRep', 12, @isnumeric); %repeat number
ip.addParameter('xyzLim', [-15 15;-15 15;-15 15], @isnumeric); %figure axis range
ip.addParameter('viewAng', [-50 10], @isnumeric); %figure view angle
ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/cytoskeleton';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(dirMod);
%--------------------------------------------------------------------------
nStep=ip.Results.nStep;
nRep=ip.Results.nRep;
xyzLim=ip.Results.xyzLim;
viewAng=ip.Results.viewAng;
%--------------------------------------------------------------------------
rec=ComRecord([dirRoot filesep 'data'],...
              [dirRoot filesep 'fig'],...
              [dirRoot filesep 'init'],'nRep',nRep,'deleteFiles',true);
nPar=numel(rec.dir_all);
%--------------------------------------------------------------------------
%Compute cytoskeleton
%--------------------------------------------------------------------------
%%
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(10,300)); %baseline length 1000;
s=ModSubstrate(0.01,4);
m=ModMembrane(true,3,0,'unit',u);
a=ModActin('unit',u);
nLas17=200;
rdnAng=(rand(nLas17,1)-0.5)*2*pi*2;
coordAssigned=[zeros(nLas17,1)-0.2 (1+(rand(nLas17,1)-0.5)*2*0.5).*[cos(rdnAng) sin(rdnAng)]];
coordAssigned=coordAssigned(abs(coordAssigned(:,3))<2,:);
p=ModFreeParticle(0,'unit',u,'coordAssigned',coordAssigned,'name','Las17'); %as wasp
p.prop = {'Particle','needFollow'};
%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{m,s},{a,m}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],m,s,a,p);
    %adjust parmeters for this particular application
    %------------------------------------------
    M.mod{M.i_mod.ModMembrane}.pm.k_V=4; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_A=8; %4
    M.mod{M.i_mod.ModMembrane}.pm.k_a=0; %12
    M.mod{M.i_mod.ModMembrane}.pm.P=0.; %0
    M.mod{M.i_mod.ModMembrane}.pm.dt=0.0002;
    M.mod{M.i_mod.ModMembrane}.pm.k_c=5; %5
    M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.2; %0.2
    M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.001;
    M.mod{M.i_mod.ModMembrane}.pm.kBT=0.0; 
    M.mod{M.i_mod.ModMembrane}.pm.V0=1223; 
    M.mod{M.i_mod.ModMembrane}.pm.A0=553; 
    M.mod{M.i_mod.ModMembrane}.pm.nAVmean=4;
    M.TypForce.pm.k_ModMembrane_ModSubstrate=4;
    %------------------------------------------
%--------------------------------------------------------------------------
%%
M.mod{M.i_mod.ModActin}.pm.konG=5;
M.mod{M.i_mod.ModActin}.pm.kbr=1;
M.mod{M.i_mod.ModActin}.pm.ksev=0.05;
M.mod{M.i_mod.ModActin}.pm.kcap=0.5;
nt=5000;
for it=1:nt
    fprintf('finishded  %d till %d\n',it,nt);
    Msave=M;
    M.mod{M.i_mod.ModFreeParticle}.var.coord(:,1)=min(M.mod{M.i_mod.ModActin}.var.coord(:,1));
    konG=obj.pm.konG*ones(M.mod{M.i_mod.ModActin}.var.Nbr,1);
    konG(M.mod{M.i_mod.ModActin}.var.coord(:,1)<M.mod{M.i_mod.ModFreeParticle}.var.coord(1,1))=0;
     M.mod{M.i_mod.ModActin} = Polymerization(M.mod{M.i_mod.ModActin},'konG',konG);
%      M.mod{M.i_mod.ModActin} = Depolymerization(M.mod{M.i_mod.ModActin});
     M.mod{M.i_mod.ModActin} = Branching(M.mod{M.i_mod.ModActin},M.mod{M.i_mod.ModFreeParticle},Mesh);
%      M.mod{M.i_mod.ModActin} = Severing(M.mod{M.i_mod.ModActin});
     M.mod{M.i_mod.ModActin} = Capping(M.mod{M.i_mod.ModActin});
end
%%
f=figure;
plot(M.mod{M.i_mod.ModActin},'f',f);rlim=[-2 2];
xlim([min(M.mod{M.i_mod.ModActin}.var.coord(:,1)) max(M.mod{M.i_mod.ModActin}.var.coord(:,1))]);
ylim([min(M.mod{M.i_mod.ModActin}.var.coord(:,2)) max(M.mod{M.i_mod.ModActin}.var.coord(:,2))]);
zlim([min(M.mod{M.i_mod.ModActin}.var.coord(:,3)) max(M.mod{M.i_mod.ModActin}.var.coord(:,3))]);
plot(M.mod{M.i_mod.ModFreeParticle},'f',f);
%%
% plot(M.mod{M.i_mod.ModFreeParticle},'f',f);
% plot(M.mod{M.i_mod.ModMembrane},'f',f);
[M.mod{M.i_mod.ModActin}] = Polymerization(M.mod{M.i_mod.ModActin});
f=figure;
plot(M.mod{M.i_mod.ModActin},'f',f);
%%
[M.mod{M.i_mod.ModActin}] = Branching(M.mod{M.i_mod.ModActin});
f=figure;
plot(M.mod{M.i_mod.ModActin},'f',f);
%%
%--------------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.var.idMesh = ...
ComMath.getMeshID(M.mod{M.i_mod.ModMembrane}.var.coord,M.Mesh.coord, M.Mesh.range,M.Mesh.d);
%--------------------------------------------------------------------------
%%
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
    dirPar{iPar}=rec.dir_all{iPar};
end
%--------------------------------------------------------------------------
%dynamics
fprintf('computing Endocytosis morphology...\n');
parfor iPar=1:nPar
% for iPar=1:1    
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModMembrane','ModMembrane_ModSubstrate'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    %stop simulation at every Tcutoff to save results
    Tcutoff=M.mod{M.i_mod.ModMembrane}.pm.dt;
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'sMax',10);
        i_loop_done=i_loop_done+CutOff.nS;
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
        %--------------------------------------------------------------------------
    end
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('finishded at i_par %d \n',iPar);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end