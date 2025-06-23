function [m] = densTran(m,k,varargin)
%--------------------------------------------------------------------------
        % densTran performs the diffusion of lipid between vertices in
        % edges for @ModMembrane
        % input: 
        % m - a @ModMembrane object
        % k - diffusion coefficient
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('k', @(x) isnumeric(x));
ip.addParameter('mex_avail', true, @islogical);
ip.addParameter('nt', 10, @isnumeric);
ip.addParameter('id_on_edg', [], @isnumeric);
ip.parse(m,k,varargin{:});
%----------------------------------------------------------------------------------------
nt=ip.Results.nt;
id_on_edg=ip.Results.id_on_edg;
if isempty(id_on_edg)
    id_on_edg=(1:m.var.n_edg)';
end
n_on_edg=numel(id_on_edg);
%----------------------------------------------------------------------------------------
%%
% m=mod.mod{1};
% m.var.dens=zeros(m.var.n_coord,1);
% % m.var.dens((m.var.coord(:,2)>8) & (m.var.coord(:,1)<-7))=10;
% xb=1; 
% r=sqrt(sum(m.var.coord.^2,2));
% idMid=(abs(m.var.coord(:,1))<xb) & (abs(m.var.coord(:,2)./r)>0.2);
% % idMid=(m.var.coord(:,2)<7)&(abs(m.var.coord(:,1))<7);
% idAll=(1:m.var.n_coord);
%     idTem=idAll(idMid);
%     idMid=randsample(idTem,80)';
%     id_on_edg=(1:m.var.n_edg)';
%     id_on_edg(sum(ismember(m.var.edge_all,idMid),2)>0)=[];
%     id_off_edg=(1:m.var.n_edg)';
%     id_off_edg(sum(ismember(m.var.edge_all,idMid),2)==0)=[];
%     
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(1,3,1)
% plot(m,'f',fig,'LineStyle','-','col',m.var.dens);
% xlim([-20 20]);ylim([-20 20]);zlim([-20 20]);view([0 90]);hold on
% for i=1:numel(id_off_edg)
%     plot3([m.var.coord(m.var.edge_all(id_off_edg(i),1),1),m.var.coord(m.var.edge_all(id_off_edg(i),2),1)],...
%           [m.var.coord(m.var.edge_all(id_off_edg(i),1),2),m.var.coord(m.var.edge_all(id_off_edg(i),2),2)],...
%           [m.var.coord(m.var.edge_all(id_off_edg(i),1),3),m.var.coord(m.var.edge_all(id_off_edg(i),2),3)],'linewidth',2);hold on;
% end
% subplot(1,3,3)
% histogram(m.var.dens)
% hold on

if ip.Results.mex_avail==false
for i=1:nt
   dDens1=k*m.var.dens(m.var.edge_all(:,1))./m.var.n_node(m.var.edge_all(:,1))';
   dDens2=k*m.var.dens(m.var.edge_all(:,2))./m.var.n_node(m.var.edge_all(:,2))';
   for iedg=1:n_on_edg
       j=id_on_edg(iedg);
       m.var.dens(m.var.edge_all(j,1))=m.var.dens(m.var.edge_all(j,1))-dDens1(j)+dDens2(j);
       m.var.dens(m.var.edge_all(j,2))=m.var.dens(m.var.edge_all(j,2))-dDens2(j)+dDens1(j);
   end
end
else
    pmc=[k,nt];
    m.var.dens=ModMembrane.densTranMex(m.var.dens,m.var.edge_all,m.var.n_node,pmc,id_on_edg);
end

% subplot(1,3,2)
% plot(m,'f',fig,'LineStyle','-','col',m.var.dens);
% xlim([-20 20]);ylim([-20 20]);zlim([-20 20]);view([0 90]);hold on;
% % for i=1:numel(id_off_edg)
% %     plot3([m.var.coord(m.var.edge_all(id_off_edg(i),1),1),m.var.coord(m.var.edge_all(id_off_edg(i),2),1)],...
% %           [m.var.coord(m.var.edge_all(id_off_edg(i),1),2),m.var.coord(m.var.edge_all(id_off_edg(i),2),2)],...
% %           [m.var.coord(m.var.edge_all(id_off_edg(i),1),3),m.var.coord(m.var.edge_all(id_off_edg(i),2),3)],'linewidth',2);hold on;
% % end
% subplot(1,3,3)
% histogram(m.var.dens)
% hold on