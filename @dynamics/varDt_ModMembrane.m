function [dtMin] = varDt_ModMembrane(dyn,f,fTotonMembrane,m,varargin)
%--------------------------------------------------------------------------
        % varDt_ModMembrane computes the adaptive time step needed by
        % @ModMembrane for dynamics
        % input: 
        % dyn - a @dynamics object
        % f - a @TypForce object
        % M - a given @Model object
        % fTotonMembrane - total force on @ModMembrane
        % m - a @ModMembrane object
        % optional:
        % see variable arguments
        %   See also TimeEval
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dyn', @(x) isa(x,'dynamics'));
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('fTotonMembrane', @(x) isnumeric(x)); 
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.parse(dyn,f,fTotonMembrane,m,varargin{:});
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%========================================================================================================================== compute dt
%%
dt=[];
%--------------------------------------------------------------------------------------------------------membrane
    d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
    r = sqrt(sum(d.^2,2));
    i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
    i = floor(r/m.pm.dr+0.5)-i_shift;
    r=floor(r/m.pm.dr+0.5)*m.pm.dr;
    c_1=(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-2)+f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-1))).^2-r.^2;
    c_2=(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i))  +f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)+1))).^2-r.^2;
      
    f_mod=fTotonMembrane*m.pm.mu;
    b = 2*(m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:)-m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:))...
        .*(f_mod(m.var.edge_all(m.var.id_on_edg,2),:)-f_mod(m.var.edge_all(m.var.id_on_edg,1),:));
    b=sum(b,2);
    dt_det=[c_1 c_2]./b;
%     dt_rand=0.5*[f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)+1)-f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)),...
%              f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-1)-f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-2)];
%     dt_rand=0.5*[f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)+1)-f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-2)];    
%     dt_rand=dt_rand.^2;
%     dt_rand=dt_rand/(2*m.pm.mu*m.pm.kBT);
     id_keep=dt_det>0;
%     dt_final=[sum(dt_det.*id_keep,2),sum(dt_rand.*id_keep,2)]; 
%     dt_final=[sum(dt_det.*id_keep,2),dt_rand]; 
    %%
%     dt_final=min(dt_final,[],2);
    dt_final=sum(dt_det.*id_keep,2);
    
    dt=[dt;min(dt_final)];
%    [min(sum(dt_det.*id_keep,2)), min(sum(dt_rand.*id_keep,2))]
%--------------------------------------------------------------------------------------------------------
dtMin=min(dt);