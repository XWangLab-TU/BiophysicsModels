function [m] = locDynLattice(m,lc,M,varargin)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('lc', @(x) isa(x,'ComLattice'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addParameter('k', 2, @isnumeric);
ip.addParameter('n', 1000, @isnumeric);
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('mex_avail', true, @islogical);
ip.parse(m,lc,M,varargin{:});
%----------------------------------------------------------------------------------------
iM=M.i_mod.ModMembrane;
k=ip.Results.k;
n=ip.Results.n;
%----------------------------------------------------------------------------------------
if ip.Results.mex_avail==false
idAlt=zeros(6,1);
idAlt(1)=lc.geometry.Dx;
idAlt(2)=-lc.geometry.Dx;
idAlt(3)=lc.geometry.Dy;
idAlt(4)=-lc.geometry.Dy;
idAlt(5)=lc.geometry.Dz;
idAlt(6)=-lc.geometry.Dz;
%----------------------------------------------------------------------------------------
for i=1:n
iVer=randsample(m.var.id_on_coord,1);
meshID1=lc.component{iM}.meshID(iVer);
meshID2=lc.component{iM}.meshID(m.var.j_T(iVer,1:m.var.n_node(iVer)));
[~,D]=distance(lc,meshID1,meshID2);
V=k*sum((D-m.pm.l0).^2);
%--------------------------------------------------------------------------
Valt=zeros(6,1);
for iAlt=1:6
meshID1=meshID1+idAlt(iAlt);
[~,D]=distance(lc,meshID1,meshID2);
Valt(iAlt)=k*sum((D-m.pm.l0).^2);
meshID1=meshID1-idAlt(iAlt);
end
%----------------------------------------------------------------------------------------
[Vmin,IDmin]=min(Valt);
if Vmin<V
    meshID1=meshID1+idAlt(IDmin);
    lc.component{iM}.meshID(iVer)=meshID1;
end
end
%----------------------------------------------------------------------------------------
m.var.coord=meshToCoord(lc,lc.component{iM}.meshID);
else
%======================================================================================== mex
pmc=zeros(10,1);
    pmc(1) = 100000;
    pmc(2) = m.pm.dt;
    pmc(3) = m.pm.P;
    pmc(4) = m.pm.k_c*0.;
    pmc(5) = m.pm.k_e;
    pmc(6) = m.pm.dr;
    pmc(7) = 0.; %m.pm.Vvol.k_V;
    pmc(8) = 0.; %m.pm.Vsurf.k_A;
    pmc(9) = m.pm.Vvol.V0*0;
    pmc(10) = m.pm.Vsurf.A0*0;
    pmc(11) = m.pm.nAVmean;
    pmc(12) = m.pm.k_a;   
    pmc(13) = m.pm.l0;   
    pmc(14) = m.pm.kBT;   
    pmc(15) = 1;  %-1:GLOBAL relax, to judge if remesh needed; 1: local relax
    j_T = m.var.j_T; j_T(isnan(j_T)) = 0;
    pmVedg=zeros(9,1);
    pmVedg(1)=m.pm.Vedg.V_0;
    pmVedg(2)=m.pm.Vedg.r_1;
    pmVedg(3)=m.pm.Vedg.r_2;
    pmVedg(4)=m.pm.Vedg.rb_1;
    pmVedg(5)=m.pm.Vedg.rb_2;
    pmVedg(6)=m.pm.Vedg.k_w;
    pmVedg(7)=m.pm.Vedg.e_b;
    pmVedg(8)=m.pm.Vedg.e_w;
    pmVedg(9)=m.pm.Vedg.k_b; 
[m.var.coord]=ModMembrane.locDynLatticeMex...
       (m.var.coord,...
        pmc,...
        m.var.edge_all,...
        m.var.face_unq,...
        j_T,...
        m.var.id_on_coord,...
        m.var.n_node',...
        m.var.T_s',...
        m.var.T_e',...
        m.var.dens,...
        pmVedg,...
        iVerNew);
%  plot(m,'FaceAlpha', 1, 'LineStyle','--');   
%========================================================================================
end
end