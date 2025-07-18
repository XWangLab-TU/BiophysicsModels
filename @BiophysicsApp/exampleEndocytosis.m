function [] = exampleEndocytosis(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % exampleEndocytosis performs the example of forming endocytosis
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, membraneFussion, exampleLamellipodia
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
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
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Endocytosis';
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
%Compute Endocytosis
%--------------------------------------------------------------------------
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); %baseline length 1000;
s=ModSubstrate(0.01,4);
m=ModMembrane(true,3,0,'unit',u);
%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{m,s}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],m,s);
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