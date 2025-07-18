function [var] = setVar(obj,update,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModClathrin'));
ip.addRequired('update', @(x) islogical(x));
ip.addParameter('n', 1, @isnumeric);
ip.addParameter('coord_org', [], @isnumeric);
ip.parse(obj,update,varargin{:});
%----------------------------------------------------------------------------------------
n=ip.Results.n;
coord_org=ip.Results.coord_org;
if (update==false) && ((isempty(n))||(isempty(coord_org)))
    error('need both n and coord_org if not updating')
end
%----------------------------------------------------------------------------------------
if update==false
    var=struct('coord',zeros(n,3),'ang',[],'n_coord',n,'coord_org',coord_org,'a',[],'ang_a',[],'O',[],'E',[],'V',[],...
        'connect',zeros(3,2,n),'ctrl_pt_org',zeros(6,3),'bound',true(n,3),'idMesh',[]);
    %var.bound(1,2)=false;
    %var.coord(:,1)=var.coord(:,1)-2;
    %----------------------------------------------------------------------
%     var.ctrl_pt_org(1,2)=0;
%     z_tem=1;
%     x_tem=(0:0.001:10)';
%     r_tem=[x_tem,zeros(size(x_tem)),ones(size(x_tem))*z_tem];
%     d_tem=abs(sqrt(sum((r_tem).^2,2))-sqrt(sum((r_tem-var.coord_org(8,1:3)).^2,2)));
%     [~,id_tem]=min(d_tem);
%     var.ctrl_pt_org(1,1)=x_tem(id_tem);
%     var.ctrl_pt_org(1,3)=z_tem;
%     var.ctrl_pt_org(2,:)=obj.rotation(var.ctrl_pt_org(1,:),0,2*pi/3);
%     var.ctrl_pt_org(3,:)=obj.rotation(var.ctrl_pt_org(2,:),0,2*pi/3);
%     var.ctrl_pt_org(4,:)=var.coord_org(4,1:3);
%     var.ctrl_pt_org(5,:)=var.coord_org(12,1:3);
%     var.ctrl_pt_org(6,:)=var.coord_org(20,1:3);    
% %     var.ctrl_pt_org(4,1)=var.coord_org(8,1)-var.ctrl_pt_org(1,1);
% %     var.ctrl_pt_org(4,2)=0;
% %     var.ctrl_pt_org(4,3)=var.coord_org(8,3)-var.ctrl_pt_org(1,3);
% %     
% %     var.ctrl_pt_org(5,:)=obj.rotation(var.ctrl_pt_org(4,:),0,2*pi/3);
% %     var.ctrl_pt_org(6,:)=obj.rotation(var.ctrl_pt_org(5,:),0,2*pi/3);
    %----------------------------------------------------------------------
    %%
    var.ctrl_pt_org(1,:)=var.coord_org(4,1:3)+[0 1 0];
    var.ctrl_pt_org(2,:)=obj.rotation(var.ctrl_pt_org(1,:),0,2*pi/3);
    var.ctrl_pt_org(3,:)=obj.rotation(var.ctrl_pt_org(2,:),0,2*pi/3);
    var.ctrl_pt_org(4,:)=var.coord_org(4,1:3)-[0 1 0];
    var.ctrl_pt_org(5,:)=obj.rotation(var.ctrl_pt_org(4,:),0,2*pi/3);
    var.ctrl_pt_org(6,:)=obj.rotation(var.ctrl_pt_org(5,:),0,2*pi/3);
    %%
    %----------------------------------------------------------------------
    
    var.ang=struct('phi',2*pi*rand(n,1),'theta',2*pi*rand(n,1),'psi',2*pi*rand(n,1));
    %var.ang_a=struct('Phi',pi*rand(n,1),'theta_a',pi*rand(n,1),'phi_a',2*pi*rand(n,1));
    var.ang_a=struct('Phi',0.0001,'theta_a',pi*rand(n,1),'phi_a',2*pi*rand(n,1));
    var.a=obj.get_a_from_ang(var.ang_a.Phi,var.ang_a.theta_a,var.ang_a.phi_a);
    var.O=obj.Omega(var.a,var.ang_a.Phi);
    %var.coord=[rand(var.n_coord,1),rand(var.n_coord,1),rand(var.n_coord,1)]*obj.pm.l0;
    %var.coord=[zeros(var.n_coord,1),zeros(var.n_coord,1),ones(var.n_coord,1)*(-min(var.coord_org(:,3)))];
    var.E=obj.getE(var.a,var.ang_a.Phi);
else
    var=obj.var;
    var.O=obj.Omega(var.a,var.ang_a.Phi);
    var.E=obj.getE(var.a,var.ang_a.Phi);
end

