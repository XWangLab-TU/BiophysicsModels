function [lc,m] = locDynLattice(m,lc,iM,varargin)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('lc', @(x) isa(x,'ComLattice'));
ip.addRequired('iM', @(x) isnumeric(x)); %which cell in lc is the membrane 
ip.addParameter('k', 1, @isnumeric);
ip.addParameter('n', 1000, @isnumeric);
ip.addParameter('mex_avail', false, @islogical);
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(m,lc,iM,varargin{:});
%----------------------------------------------------------------------------------------
k=ip.Results.k;
n=ip.Results.n;
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
end