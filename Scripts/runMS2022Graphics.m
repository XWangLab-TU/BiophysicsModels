%==========================================================================
% This script adjust the graphic output for the 2022 manuscript. 
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
addpath(genpath(dirMod));
%% ==========================================================================
% Graphics
%==========================================================================
% Adjust figures and make a movie in a given directory
%% 1.Red blood cell
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig3_RBC/fig/rep01';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-10 10;-10 10;-10 10];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(5).LineWidth=0.5; A.Children(5).EdgeColor=[0 1 1];
        A.XTick=(-10:5:10); A.YTick=(-10:5:10); A.ZTick=(-10:5:10);
        view([16 47]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.005,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the last figure larger
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
%% 2.Vesicle Fussion
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig3_fussion/fig/rep01';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-5 11;-8 8;-8 8];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(5).LineWidth=0.5; A.Children(5).EdgeColor=[0 1 1];
        A.XTick=(-5:5:10); A.YTick=(-5:5:5); A.ZTick=(-5:5:5);
        view([15 12]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.005,'n_file',40,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
%make the 1st, 2nd and last figure larger
myDir=dir(dirFig);
for i=3:4
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
for i=7:8
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
nDir=numel(myDir);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
%% 3.Filopodia
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Filopodia/fig/rep03';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-11 17;-12 16;-14 14];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-10:5:15); A.YTick=(-10:5:15); A.ZTick=(-10:10:10);
        view([80 62]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.005,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the last figure larger
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%% 4.Lamellipodia
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Lamellipodia/fig/rep08';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-10 16;-13 13;-13 13];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-10:5:15); A.YTick=(-10:5:15); A.ZTick=(-10:10:10);
        view([10 65]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.005,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the last figure larger
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%% 5.Endocytosis
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Endocytosis/fig/rep01';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-9 9;-9 9;-9 9];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-5:5:5); A.YTick=(-5:5:5); A.ZTick=(-5:5:5);
        view([-56 2]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.005,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the last figure larger
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%% 6.Thether pulling
%normal tether
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/normalTether/fig/rep07';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-15 15;-10 20;-15 15];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-15:5:15); A.YTick=(-5:5:15); A.ZTick=(-15:5:15);
        view([0 90]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.0005,'print_info','%.4f','view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%extra membrane tether
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/extraMemTether/fig/rep10';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-15 15;-10 20;-15 15];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-15:5:15); A.YTick=(-5:5:15); A.ZTick=(-15:5:15);
        view([0 90]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.0005,'print_info','%.4f','view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%tether pulling: normal membrane
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/fig/p04_rep1';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-40 15;-10 45;-15 15];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-35:10:15); A.YTick=(-10:10:40); A.ZTick=(-15:10:15);
        view([0 90]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.0005,'print_info','%.4f','view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the 1st and last figure larger
mkdir([dirFig filesep 'Large']);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
for i=5:6
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%tether pulling: diffusion barrier
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/fig/p10_rep6';
myDir=dir(dirFig);
nDir=numel(myDir);
dirData='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/data/p10_rep6';
myDirData=dir(dirData);
nDirData=numel(myDirData);
%adjust every figure nicely
xyzLim=[-40 15;-10 45;-15 15];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-35:10:15); A.YTick=(-10:10:40); A.ZTick=(-15:10:15);
        view([0 90]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.0005,'print_info','%.4f','view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the 1st and last figure larger, also show diffusion barrier
mkdir([dirFig filesep 'Large']);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        S=load([myDirData(end).folder filesep myDirData(end).name]);
        M=S.M;
        xb=M.mod{M.i_mod.ModMembrane}.pm.xbDiff;
        r_thre=M.mod{M.i_mod.ModMembrane}.pm.rDiff;
        r=sqrt(sum(M.mod{M.i_mod.ModMembrane}.var.coord.^2,2));
        idMid=(abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,1))<xb) & (abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,2)./r)>r_thre);
        idAll=(1:M.mod{M.i_mod.ModMembrane}.var.n_coord);
        idMid=idAll(idMid);
        id_on_edg=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
        id_barrier=id_on_edg(sum(ismember(M.mod{M.i_mod.ModMembrane}.var.edge_all,idMid),2)>0);
        n_barrier=numel(id_barrier);
        for j=1:n_barrier
            hold on;
            edg=M.mod{M.i_mod.ModMembrane}.var.edge_all(id_barrier(j),:);
            plot3([M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),1),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),1)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),2),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),2)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),3),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),3)],...
                'linewidth',2,'color',[1 1 0]);
        end
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
for i=3:6
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        S=load([myDirData(3).folder filesep myDirData(3).name]);
        M=S.M;
        xb=M.mod{M.i_mod.ModMembrane}.pm.xbDiff;
        r_thre=M.mod{M.i_mod.ModMembrane}.pm.rDiff;
        r=sqrt(sum(M.mod{M.i_mod.ModMembrane}.var.coord.^2,2));
        idMid=(abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,1))<xb) & (abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,2)./r)>r_thre);
        idAll=(1:M.mod{M.i_mod.ModMembrane}.var.n_coord);
        idMid=idAll(idMid);
        id_on_edg=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
        id_barrier=id_on_edg(sum(ismember(M.mod{M.i_mod.ModMembrane}.var.edge_all,idMid),2)>0);
        n_barrier=numel(id_barrier);
        for j=1:n_barrier
            hold on;
            edg=M.mod{M.i_mod.ModMembrane}.var.edge_all(id_barrier(j),:);
            plot3([M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),1),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),1)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),2),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),2)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),3),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),3)],...
                'linewidth',2,'color',[1 1 0]);
        end
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%tether pulling: diffusion barrier+extra membrane
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/fig/p16_rep1';
myDir=dir(dirFig);
nDir=numel(myDir);
dirData='/endosome/work/bioinformatics/s171152/data/membrane/Fig5/data/p16_rep1';
myDirData=dir(dirData);
nDirData=numel(myDirData);
%adjust every figure nicely
xyzLim=[-40 15;-10 45;-15 15];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(1).SizeData=5; A.Children(6).LineWidth=0.5; A.Children(6).EdgeColor=[0 1 1];
        A.XTick=(-35:10:15); A.YTick=(-10:10:40); A.ZTick=(-15:10:15);
        view([0 90]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 5 5]);
        close(gcf);
    end
end
ComRecord.makeStackMov(dirFig,'DelayTime',0.1,'dt',0.0005,'print_info','%.4f','view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
%make the 1st and last figure larger, also show diffusion barrier
mkdir([dirFig filesep 'Large']);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        S=load([myDirData(end).folder filesep myDirData(end).name]);
        M=S.M;
        xb=M.mod{M.i_mod.ModMembrane}.pm.xbDiff;
        r_thre=M.mod{M.i_mod.ModMembrane}.pm.rDiff;
        r=sqrt(sum(M.mod{M.i_mod.ModMembrane}.var.coord.^2,2));
        idMid=(abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,1))<xb) & (abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,2)./r)>r_thre);
        idAll=(1:M.mod{M.i_mod.ModMembrane}.var.n_coord);
        idMid=idAll(idMid);
        id_on_edg=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
        id_barrier=id_on_edg(sum(ismember(M.mod{M.i_mod.ModMembrane}.var.edge_all,idMid),2)>0);
        n_barrier=numel(id_barrier);
        for j=1:n_barrier
            hold on;
            edg=M.mod{M.i_mod.ModMembrane}.var.edge_all(id_barrier(j),:);
            plot3([M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),1),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),1)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),2),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),2)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),3),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),3)],...
                'linewidth',2,'color',[1 1 0]);
        end
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
for i=3:6
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        S=load([myDirData(3).folder filesep myDirData(3).name]);
        M=S.M;
        xb=M.mod{M.i_mod.ModMembrane}.pm.xbDiff;
        r_thre=M.mod{M.i_mod.ModMembrane}.pm.rDiff;
        r=sqrt(sum(M.mod{M.i_mod.ModMembrane}.var.coord.^2,2));
        idMid=(abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,1))<xb) & (abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,2)./r)>r_thre);
        idAll=(1:M.mod{M.i_mod.ModMembrane}.var.n_coord);
        idMid=idAll(idMid);
        id_on_edg=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
        id_barrier=id_on_edg(sum(ismember(M.mod{M.i_mod.ModMembrane}.var.edge_all,idMid),2)>0);
        n_barrier=numel(id_barrier);
        for j=1:n_barrier
            hold on;
            edg=M.mod{M.i_mod.ModMembrane}.var.edge_all(id_barrier(j),:);
            plot3([M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),1),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),1)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),2),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),2)],...
                  [M.mod{M.i_mod.ModMembrane}.var.coord(edg(1),3),M.mod{M.i_mod.ModMembrane}.var.coord(edg(2),3)],...
                'linewidth',2,'color',[1 1 0]);
        end
        ComRecord.savePlot([dirFig filesep 'Large'],[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 8]);
    end
end
close(gcf);
%adjust analysis figures
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/Fig5';
myDir=dir(dirFig);
nDir=numel(myDir);
for i=3:nDir
    if numel(myDir(i).name) > 3
    if strcmp(myDir(i).name(end-3:end),'.fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig ,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 8 6]);
        close(gcf);
    end
    end
end
%% 7.Remeshing test
dirFig='/endosome/work/bioinformatics/s171152/data/membrane/FigS3_SupTest/CS/fig/rep02';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-5 11;-8 8;-8 8];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(end).FaceAlpha=0.25;
%         A.Children(7).LineWidth=0.1; A.Children(7).EdgeColor=[0 1 1];
        A.Children(end).LineStyle='none';
        A.XTick=(-5:5:10); A.YTick=(-5:5:5); A.ZTick=(-5:5:5);
        view([15 12]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.2,'dt',0.005,'n_file',40,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
%make the 1st, 2nd and last figure larger
myDir=dir(dirFig);
for i=3:4
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
for i=7:8
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
nDir=numel(myDir);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);

dirFig='/endosome/work/bioinformatics/s171152/data/membrane/FigS3_SupTest/ORG/fig/rep01';
myDir=dir(dirFig);
nDir=numel(myDir);
%adjust every figure nicely
xyzLim=[-5 11;-8 8;-8 8];
for i=3:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        A=gca;
        A.Children(end).FaceAlpha=0.25;
%         A.Children(7).LineWidth=0.1; A.Children(7).EdgeColor=[0 1 1];
        A.Children(end).LineStyle='none';
        A.XTick=(-5:5:10); A.YTick=(-5:5:5); A.ZTick=(-5:5:5);
        view([15 12]);
        material shiny;
        ComRecord.savePlot(dirFig,myDir(i).name(1:end-4),'figH',gcf,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
        close(gcf);
    end
end
%make a gif movie
ComRecord.makeStackMov(dirFig,'DelayTime',0.2,'dt',0.005,'n_file',40,'view_ang',[],'xyzLim',xyzLim,'PaperPosition', [0 0 7.5 5],'OuterPosition', [0 0 7.5 5]);
%make the 1st, 2nd and last figure larger
myDir=dir(dirFig);
for i=3:4
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
for i=7:8
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);
myDir=dir(dirFig);
nDir=numel(myDir);
for i=nDir-1:nDir
    if strcmp(myDir(i).name(end-2:end),'fig')
        open([myDir(i).folder filesep myDir(i).name]);
        ComRecord.savePlot(dirFig,[myDir(i).name(1:end-4) 'LARGE'],'figH',gcf,'PaperPosition', [0 0 8 5.33],'OuterPosition',[0 0 8 5.33]);
    end
end
close(gcf);