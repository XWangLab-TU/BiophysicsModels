%==========================================================================
% This script includes all the examples studied in the 2022 manuscript. 
% Applications are separated into sections by double %% signs 
% and can be run individually with 'ctrl'+'enter'
% see the documentation of individual application for more detail
%   See also runMS2022Graphics
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/25
%==========================================================================
%% add library
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels'; 
% dirMod='/home2/s171152/codes/matlab/mine/git/module/module';
addpath(genpath(dirMod));
%% ========================================================================
% Computation
%==========================================================================
%% Red blood cell
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig3_RBC';
tic
BiophysicsApp.redBloodCell(dirRoot,dirMod,'nRep', 12); %repeat 12 times
toc
%% Fusion
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig3_fussion';
tic
BiophysicsApp.membraneFussion(dirRoot,dirMod,'nRep', 12);
toc
%% Filopodia
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Filopodia';
tic
BiophysicsApp.exampleFilopodia(dirRoot,dirMod,'nRep', 12);
toc
%% Lamellipodia
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Lamellipodia';
tic
BiophysicsApp.exampleLamellipodia(dirRoot,dirMod,'nRep', 12);
toc
%% Endocytosis
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Endocytosis';
tic
BiophysicsApp.exampleEndocytosis(dirRoot,dirMod,'nRep', 12);
toc
%% Tether pulling
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig5';
BiophysicsApp.tetherPull(dirRoot,dirMod,'kDiffAlt', [0.025,0.05,0.1,0.15,0.2,0.25],'reloadInit',true,'idInitNormal',6,'idInitExtraMem',10,...
                        'reloadInitExtraMem',true,'reloadComputed',false,'reloadAnalysis',false,'r_threAlt',0.1);
%% TestRemeshGaussianCurv
dirRoot= '/endosome/work/bioinformatics/s171152/data/membrane/FigS3_SupTest';
dirRoot1 = '/endosome/work/bioinformatics/s171152/data/membrane/FigS3_SupTest/CS';
dataStored=BiophysicsApp.membraneFussion(dirRoot1,dirMod,'nRep', 12,'CSremesh',true);
cd(dirRoot1);
save(dataStored,'temp.mat');
%% Test Remesh
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/FigS3_SupTest';
BiophysicsApp.membraneFussion(dirRoot,dirMod,'nRep', 12,'CSremesh',true,'nStep', 20000,'CSremeshAnalysisOnly',true);
%% Test GaussianCurv
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/FigS3B_SupTest';
BiophysicsApp.membraneFussion(dirRoot,dirMod,'nRep', 12,'nStep', 10000,'GaussianCurv',true);
%% Memory Test
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/FigS3C_SupTest';
BiophysicsApp.memoryTest(dirRoot,dirMod,'nRep', 12,'nStep', 10000,'N',3);