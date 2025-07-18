nCLmax=25;
nParaX=2;
nParaY=10;
nPara=nParaX*nParaY;
nRep=1;
nPar=12;
pathName='M2';%'CLadding';
dirFolder=cell(nCLmax,1);
temFormat3=['%0.' num2str(3) 'd'];
for iCL=1:nCLmax
    dirFolder{iCL}=[dirRoot filesep pathName filesep num2str(iCL,temFormat3)]; mkdir(dirFolder{iCL});
end
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

T=nan(nCLmax,nPara);
k=nan(nPara,1);
for iCL=2:nCLmax  
   for iPara=1:nPara
       if trueState(iCL,iPara)==true
       dirTem=dir(dirState{iCL,iPara});
       S=load([dirTem(end).folder filesep dirTem(end).name]);
       M=S.M;
       c=M.mod{M.i_mod.ModClathrin};
       [T(iCL,iPara)] = getTightness(c,M);
       end
   end
end
iCL=2;
for iPara=1:nPara
    dirTem=dir(dirState{iCL,iPara});
    S=load([dirTem(end).folder filesep dirTem(end).name]);
    M=S.M;
    f=M.TypForce;
    k(iPara)=f.pm.k_ModClathrin(1);
end
%% relax plot
Trex=[];
iCL=3; iPara=1;
dirTem=dir(dirState{iCL,iPara});
Ndir=size(dirTem,1);
for iDir=3:Ndir
    S=load([dirTem(iDir).folder filesep dirTem(iDir).name]);
    M=S.M;
    c=M.mod{M.i_mod.ModClathrin};
    [Ttem] = getTightness(c,M);
    Trex=[Trex;Ttem];
end
plot(Trex,'.-'); xlabel('Time step');
%%
Tsingle=[];
iPara=2; %M1 02
for iCL=2:nCLmax
    if trueState(iCL,iPara)==true
        Tsingle=[Tsingle;T(iCL,iPara)];
    end
end
figure;plot((2:size(Tsingle,1)+1),Tsingle,'.-');
Tsingle=[];
iPara=16; %M1 16
for iCL=2:nCLmax
    if trueState(iCL,iPara)==true
        Tsingle=[Tsingle;T(iCL,iPara)];
    end
end
figure;plot((2:size(Tsingle,1)+1),Tsingle,'.-');
%%
% iCL=25; iPara=2;
% dirTem=dir(dirState{iCL,iPara});
% S=load([dirTem(end).folder filesep dirTem(end).name]);
% xyzLimPlot=[-15 15;-10 20;-15 15];viewAng=[0 0];
% fig=ComPlot.PlotMod(S.M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
%            [S.M.i_mod.ModClathrin,S.M.i_mod.ModFreeParticle,S.M.i_mod.ModMembrane]);
%        xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:));view(viewAng);
%        axis off;
%        w=xyzLimPlot(2,2)-xyzLimPlot(2,1);
%        h=xyzLimPlot(3,2)-xyzLimPlot(3,1);
%        h=h/w;
%        ComRecord.savePlot(pwd, 'open','figH', fig,'PaperPosition',[0 0 figSize figSize*h],'OuterPosition',[0 0 figSize figSize*h]);
%        close(fig);
% %%
% iPara=16;
% for iCL=10:13
% dirTem=dir(dirState{iCL,iPara});
% S=load([dirTem(end).folder filesep dirTem(end).name]);
% xyzLimPlot=[-15 15;-10 20;-15 15];viewAng=[0 0];
% fig=ComPlot.PlotMod(S.M, 'xyzLim', xyzLimPlot,'viewAng',viewAng,'idModPlot',...
%            [S.M.i_mod.ModClathrin,S.M.i_mod.ModFreeParticle,S.M.i_mod.ModMembrane]);
%        xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:));view(viewAng);
%        axis off;
%        w=xyzLimPlot(2,2)-xyzLimPlot(2,1);
%        h=xyzLimPlot(3,2)-xyzLimPlot(3,1);
%        h=h/w;
%        ComRecord.savePlot(pwd, ['closed' num2str(iCL)],'figH', fig,'PaperPosition',[0 0 figSize figSize*h],'OuterPosition',[0 0 figSize figSize*h]);
%        close(fig);
% end
%%
Tplot=[];
Ttem=[];
for iPara=1:2:nPara
    Ttem=[Ttem;max(T(:,iPara))];
%     for iCL=nCLmax:-1:2 
%         if ~isnan(T(iCL,iPara))
%             Ttem=[Ttem;T(iCL,iPara)];
%             break;
%         end
%     end
end
Tplot=[Tplot,Ttem];
Ttem=[];
for iPara=2:2:nPara
    Ttem=[Ttem;max(T(:,iPara))];
%     for iCL=nCLmax:-1:2 
%         if ~isnan(T(iCL,iPara))
%             Ttem=[Ttem;T(iCL,iPara)];
%             break;
%         end
%     end
end
Tplot=[Tplot,Ttem];
plot(k(1:2:nPara),Tplot(:,1),'.-'); hold on;
plot(k(2:2:nPara),Tplot(:,2),'.-'); hold on;