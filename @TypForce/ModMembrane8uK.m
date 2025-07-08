function [m] = ModMembrane8uK(f,m, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(f,m, varargin{:});
%--------------------------------------------------------------------------------------------------------
K_tem = m.var.f.K./m.var.f.A;
K_n=sqrt(sum(K_tem.*K_tem,2));
m.var.f.u_K = K_tem./K_n;
%--------------------------------------------------------------------------
id_tem=sum(m.var.u_face.*m.var.f.u_K,2) < 0;
m.var.f.u_K(id_tem,:)=-m.var.f.u_K(id_tem,:);

%%
% fig=figure('units','normalized','outerposition',[0 0 1 1]); 
% view_theta=0; view_phi=50; x_lim=[-15 15];
% plot(m,'linestyle','-','f',fig,'FaceAlpha',1,'view_ang',[view_theta view_phi]);
% hold on;
% quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),m.var.f.u_K(:,1),m.var.f.u_K(:,2),m.var.f.u_K(:,3),'linewidth',2);
% % hold on;
% % quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),dir_abc_save(:,1),dir_abc_save(:,2),dir_abc_save(:,3),'linewidth',2);
% xlim(x_lim);ylim(x_lim);zlim(x_lim);
end
