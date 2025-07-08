function f = PlotMod(M, varargin)
%--------------------------------------------------------------------------
        % PlotMod plots all the objects in the given @Model object M
        % input: 
        % dirAll - directory structure for data and figure
        % M - a given @Model object
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('f', [], @isobject);
ip.addParameter('xyzLim', [], @isnumeric); %support multi-angle views as subplot with >1 in 3rd dimension
ip.addParameter('viewAng', [0 0], @isnumeric); %support multi-angle views as subplot with >1 rows
ip.addParameter('subplotId', [1 1], @isnumeric); %subplot dimension: [row, col] if multi-view assigned
ip.addParameter('idModPlot', [], @isnumeric);
ip.parse(M, varargin{:});
%==========================================================================
f = ip.Results.f;
if isempty(f)
    f=figure;
    figure(f); hold on;
else
    figure(f); hold on;
end
xyzLim=ip.Results.xyzLim;
if isempty(xyzLim)
    xyzLim=[-10 10;-10 10;-10 10];
end
viewAng=ip.Results.viewAng;
subplotId=ip.Results.subplotId;
nSubplot=1;
if numel(size(xyzLim))>2
    nSubplot=size(xyzLim,3);
    if size(viewAng,1)~=nSubplot
        error('xyzLim and viewAng no match');
    elseif isempty(subplotId)
        warning('subplotId assigned by default');
        subplotId=[1,nSubplot];
    end
end
idModPlot=ip.Results.idModPlot;
if isempty(idModPlot)
    idModPlot=1:M.n_mod;
end
nModPlot=numel(idModPlot);
%==========================================================================
for iSubPlot=1:nSubplot
    figure(f);
    subplot(subplotId(1),subplotId(2),iSubPlot);
for iModPlot=1:nModPlot
    i=idModPlot(iModPlot);
    if strcmp(M.name{i},'ModMembrane')
        plot(M.mod{M.i_mod.ModMembrane},'f',f,'LineStyle','-','facealpha',1,'LineWidth', 0.5,'LineCol',[0.1490    0.1490    0.1490]);
        a=gca;
        h = findobj(a,'Type','Patch');
        h.FaceColor=[0.5804 0.1608 0.1608];
    elseif strcmp(M.name{i},'ModFreeParticle')
        plot(M.mod{M.i_mod.ModFreeParticle},'f',f,'Color',[0 1 0]);
    elseif strcmp(M.name{i},'ModClathrin')
        col_tem=zeros(M.mod{M.i_mod.ModClathrin}.var.n_coord,3); col_tem(:,2)=0.;col_tem(:,2)=1;col_tem(:,3)=1;
        plot(M.mod{M.i_mod.ModClathrin},'f',f,'simple',true,'col',col_tem);
    elseif strcmp(M.name{i},'ModSubstrate')
        plot(M.mod{M.i_mod.ModSubstrate},'f',f,'Color',[0 1 0]);
    elseif strcmp(M.name{i},'ModPolymerChain')
        plot(M.mod{M.i_mod.ModPolymerChain},'f',f,'facealpha',1);
    else
        plot(M.mod{M.i_mod.(M.name{i})},'f',f);
    end
end
    if nSubplot<2
        xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
        lighting gouraud;
        view(viewAng);
    else
        xlim(xyzLim(1,:,iSubPlot));ylim(xyzLim(2,:,iSubPlot));zlim(xyzLim(3,:,iSubPlot));
        lighting gouraud;
        view(viewAng(iSubPlot,:));
    end
end


