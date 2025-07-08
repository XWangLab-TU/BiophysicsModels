function [f,F_tot] = ModRingmem(f,rm,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('rm', @(x) isa(x,'ModRingmem'));
ip.addParameter('mex', [], @isobject);
ip.parse(f,rm,varargin{:});
%----------------------------------------------------------------------------------------
%%
rm=SetVar(rm,'update',true);
[F_i,F_g,F_iIntp]=compF(rm);
F_tot=sum(F_i)+F_g;

f.int_comp=zeros(rm.var.n_coord,2);
% plot(rm.var.coordIntp(:,1),rm.var.coordIntp(:,2),'.-');hold on;plot(rm.var.coord(:,1),rm.var.coord(:,2),'o');hold on;
%%
for i=1:rm.var.n_coord-2
    %%
    rm_alt=rm;
    rm_alt.var.coord(i,1)=rm_alt.var.coord(i,1)+rm.pm.dR;
    rm_alt=SetVar(rm_alt,'update',true);
%     plot(rm_alt.var.coordIntp(:,1),rm_alt.var.coordIntp(:,2),'.-');hold on;plot(rm_alt.var.coord(:,1),rm_alt.var.coord(:,2),'o');hold on;
%     [F_i_alt,F_g_alt]=compF(rm_alt);
%     if i<=rm.var.n_coord-1
%         [F_i_alt,F_g_alt]=compF(rm_alt,'idx',i,'F_old',F_i,'F_iIntp_old',F_iIntp,'obj_old',rm);
%     else
        [F_i_alt,F_g_alt]=compF(rm_alt);
%     end
%     [S]=rm.ModRingmem8compS(rm.var.coord(i,1),rm.var.coord(i+1,1),rm.var.coord(i,2),rm.var.coord(i+1,2));
    f.int_comp(i,1)=-(sum(F_i_alt)-sum(F_i)+F_g_alt-F_g)/rm.pm.dR;
%--------------------------------------------------------------------------
    %%
    rm_alt=rm;
    rm_alt.var.coord(i,2)=rm_alt.var.coord(i,2)+rm.pm.dZ;
    rm_alt=SetVar(rm_alt,'update',true);
%     [F_i_alt,F_g_alt]=compF(rm_alt);
%     if i<=rm.var.n_coord-1
%         [F_i_alt,F_g_alt]=compF(rm_alt,'idx',i,'F_old',F_i,'F_iIntp_old',F_iIntp,'obj_old',rm);
%     else
        [F_i_alt,F_g_alt]=compF(rm_alt);
%     end
%     [S]=rm.ModRingmem8compS(rm.var.coord(i,1),rm.var.coord(i+1,1),rm.var.coord(i,2),rm.var.coord(i+1,2));
    f.int_comp(i,2)=-(sum(F_i_alt)-sum(F_i)+F_g_alt-F_g)/rm.pm.dR;
end
end