%% ========================================================================
nCLmax=25;
nParaX=2;
nParaY=10;
nPara=nParaX*nParaY;
nRep=1;
nPar=12;
dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
pathName='M2';%'CLadding';
dirFolder=cell(nCLmax,1);
temFormat3=['%0.' num2str(3) 'd'];
for iCL=1:nCLmax
    dirFolder{iCL}=[dirRoot filesep pathName filesep num2str(iCL,temFormat3)]; mkdir(dirFolder{iCL});
end
xyzLimPlot=[-15 15;-15 15;-15 15];viewAng=[0 0];
figSize=10;
%%
temFormat2=['%0.' num2str(2) 'd'];
trueState=true(nCLmax,nPara);
dirState=cell(nCLmax,nPara);
for iCL=2:nCLmax    
   for iPara=1:nPara
       dirState{iCL,iPara}=[dirFolder{iCL} filesep 'data' filesep 'p' num2str(iPara, temFormat2) '_rep1'];
       dirInfo=dir(dirState{iCL,iPara});
       nFile=numel(dirInfo);
       if nFile<10
           trueState(iCL,iPara)=false;
       end
   end
end
figPath=[dirRoot filesep pathName filesep 'fig']; mkdir(figPath);
dirFig=cell(nCLmax,nPara);
for iCL=1:nCLmax    
   for iPara=1:nPara
       dirFig{iCL,iPara}=[figPath filesep num2str(iCL,temFormat3) filesep num2str(iPara, temFormat2)];
       mkdir(dirFig{iCL,iPara});
   end
end
%%
disp('fig');
for iCL=2:nCLmax  
    disp(iCL);
   for iPara=1:nPara
       if trueState(iCL,iPara)==true
       dirTem=dir(dirState{iCL,iPara});
       for iFile=1:2
           if iFile==1
               S=load([dirTem(3).folder filesep dirTem(3).name]);
           else
               S=load([dirTem(end).folder filesep dirTem(end).name]);
           end
       M=S.M;
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
       %%
       fig=ComPlot.PlotMod(M2, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
           [M2.i_mod.ModClathrin,M2.i_mod.ModFreeParticle,M2.i_mod.ModMembrane]);
       xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:));view(viewAng);
       axis off;
       w=xyzLimPlot(2,2)-xyzLimPlot(2,1);
       h=xyzLimPlot(3,2)-xyzLimPlot(3,1);
       h=h/w;
       ComRecord.savePlot(dirFig{iCL,iPara}, num2str(iFile),'figH', fig,'PaperPosition',[0 0 figSize figSize*h],'OuterPosition',[0 0 figSize figSize*h]);
       close(fig);
       end      
       end
   end
end
clear S;
disp('gif');
gifPath=[dirRoot filesep pathName filesep 'gif']; mkdir(figPath);
dirGif=cell(nPara,1);
for iPara=1:nPara
    dirGif{iPara}=[gifPath filesep num2str(iPara, temFormat2)];
    mkdir(dirGif{iPara});
end

for iPara=1:nPara
    disp(iPara);
   for iCL=2:nCLmax 
       if trueState(iCL,iPara)==true
           copyfile([dirFig{iCL,iPara} filesep '2.fig'], [dirGif{iPara} filesep num2str(iCL, temFormat3) '.fig'], 'f');
       end
   end
   ComRecord.makeStackMov(dirGif{iPara},'DelayTime',0.2,'PaperPosition', [0,0,5,5],'OuterPosition', [0,0,5,5]);
end