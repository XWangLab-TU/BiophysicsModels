%% ========================================================================
temFormat3=['%0.' num2str(3) 'd'];
xyzLimPlot=[-15 15;-15 15;-15 15];viewAng=[0 0];
figSize=10;
%%
       dirTem=dir(pwd);
%        S=load([dirTem(end).folder filesep dirTem(end).name]);
       S=load([dirTem(3).folder filesep dirTem(3).name]);

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

       fig=ComPlot.PlotMod(M2, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
           [M2.i_mod.ModClathrin,M2.i_mod.ModFreeParticle,M2.i_mod.ModMembrane]);
       xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:));view(viewAng);
       axis off;
       w=xyzLimPlot(2,2)-xyzLimPlot(2,1);
       h=xyzLimPlot(3,2)-xyzLimPlot(3,1);
       h=h/w;
       %ComRecord.savePlot(pwd, 'switchSnapshot','figH', fig,'PaperPosition',[0 0 figSize figSize*h],'OuterPosition',[0 0 figSize figSize*h]);

clear S;