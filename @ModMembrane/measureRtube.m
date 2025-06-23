function [pCircle] = measureRtube(m,n,V0,varargin)
%--------------------------------------------------------------------------
        % measureRtube performs the measurement of radius of membrane tube
        % input: 
        % m - @ModMembrane object
        % n - unit direction of tube
        % V0 - center of measurement
        % output:
        % pCircle - [X of 2D center (from V0), Y of 2D center (from V0), radius]
        % optional:
        % see variable arguments
        %   See also MeshOri
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('n', @(x) isnumeric(x));
ip.addRequired('V0', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical); %whether to plot measured ring around membrane tube
ip.addParameter('fig', [], @isobject); %given figure to plot on top of
ip.parse(m,n,V0,varargin{:});
%----------------------------------------------------------------------------------------
plot_or_not=ip.Results.plot_or_not;
fig=ip.Results.fig;
%----------------------------------------------------------------------------------------
I=[];
pCircle=[];
idTem=sqrt(sum((m.var.coord-V0).^2,2))<5;
idAll=(1:m.var.n_coord);
idTem=idAll(idTem);
if ~isempty(idTem)
idEdg0=ismember(m.var.edge_all(:,1),idTem);
idEdg1=ismember(m.var.edge_all(:,2),idTem);
idEdg=idEdg0 | idEdg1;
idAll=(1:m.var.n_edg);
idEdg=idAll(idEdg);
nEdg=numel(idEdg);
edg=m.var.edge_all;
coord=m.var.coord;
for i=1:nEdg
    P0=coord(edg(idEdg(i),1),:);
    P1=coord(edg(idEdg(i),2),:);
    [Item,check]=ComMath.plane_line_intersect(n,V0,P0,P1);
    if check==1
        I=[I;Item];
    end
end
if ~isempty(I)
if n(1)>0
theta=0.5*pi;
phi=0.25*pi;
else
theta=0.5*pi;
phi=-0.25*pi;
end
Inew=I-mean(I);
Inew = ComMath.rotAxis(Inew',0,phi); Inew=Inew';
Inew = ComMath.rotAxis(Inew',theta,0); Inew=Inew';
if size(Inew,1)>10
pCircle = ComMath.fit2Dcircle(Inew(:,1:2));
end
if plot_or_not==true
    if isempty(fig)
        %%
fig=plot(m);
figure(fig); 
scatter3(I(:,1),I(:,2),I(:,3),40,'filled','markerfacecolor',[0 0 1]); hold on;
for i=1:nEdg
    plot3([coord(edg(idEdg(i),1),1),coord(edg(idEdg(i),2),1)],...
          [coord(edg(idEdg(i),1),2),coord(edg(idEdg(i),2),2)],...
          [coord(edg(idEdg(i),1),3),coord(edg(idEdg(i),2),3)],'-','linewidth',2,'Color',[1 0 0]); hold on;
end
figure(); 
Iold=I-mean(I);
scatter3(Iold(:,1),Iold(:,2),Iold(:,3),40,'filled','markerfacecolor',[1 0 0]); hold on;
scatter3(Inew(:,1),Inew(:,2),Inew(:,3),40,'filled','markerfacecolor',[0 0 1]); hold on;
viscircles(pCircle(1:2),pCircle(3));
xlabel('x'); ylabel('y'); zlabel('z');
    else
        figure(fig); hold on;
        scatter3(I(:,1),I(:,2),I(:,3),40,'filled','markerfacecolor',rand(3,1)); hold on;
    end
end
end
end