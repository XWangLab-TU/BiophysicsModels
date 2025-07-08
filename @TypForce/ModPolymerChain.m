function [f,Vtot] = ModPolymerChain(f,M,varargin)
%--------------------------------------------------------------------------
        % ModPolymerChain performs the computation of the forces in
        % @ModPolymerChain, including LennardJones, FENE etc
        % input: 
        % f - @TypForce
        % md - @ModPolymerChain 
        % output:
        % Vtot - total potential in @ModPolymerChain 
        % optional:
        % see variable arguments
        %   See also ModMembrane
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/18
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(f,M,varargin{:});
%----------------------------------------------------------------------------------------
md=M.mod{M.i_mod.ModPolymerChain};
%----------------------------------------------------------------------------------------
%% compute constant force, including LJ and FENE
%----------------------------------------------------------------------------------------
Mname='ModPolymerChain';
if isempty(f.int_stored.(Mname))
    %define parameters for LJ and FENE potential
    Vname='ModPolymerChain';
    Vin_pm(1)=md.pm.Vpm.epsilon;  Vin_pm(2)=md.pm.Vpm.sigma; Vin_pm(3)=md.pm.Vpm.LJorder; % for LJ potential
    Vin_pm(4)=md.pm.Vpm.kFENE;  Vin_pm(5)=md.pm.Vpm.R0;  % for FENE potential
    f=storedForce(f,Vname,Vin_pm,Mname,[md.pm.Vpm.r_1 md.pm.Vpm.dr md.pm.Vpm.r_2],'rsq_std', md.pm.Vpm.rsq_std, 'std_std',md.pm.Vpm.std_std);
end
idPair=md.var.edg;
idCoord=1:md.var.n_coord;
[f,Vtot] = constForce(f,Mname,md,idPair,idCoord);
%----------------------------------------------------------------------------------------
%% compute regular force, from worm like potential
%----------------------------------------------------------------------------------------
f.int_comp.(Mname)=zeros(md.var.n_coord,3);
Vworm=nan(md.var.n_coord,1);
VwormDir=nan(md.var.nChain,1);
Nmark=[0;cumsum(md.var.nSubunit)]; %cumsum to get the last subunit number at the end of each chain
[md] = getCosAng(md);
Vin_pm=[md.pm.Vpm.kappa,md.pm.kBT,md.pm.Vpm.Theta0]; % for worm potential
for iC=1:md.var.nChain
    iS=Nmark(iC)+1:Nmark(iC+1);
    [Vworm(iS),~] = potential.Vworm(md.var.cosAng(iS),Vin_pm);
    VwormDir(iC)=potential.VwormDir(md.var.dir(iS,:),md.pm.Vpm.kappaDir);
end
dr=md.pm.Vpm.dr;
VwormAlt=Vworm;
for iC=1:md.var.nChain
    Vtot=Vtot+sum(Vworm(Nmark(iC)+2:Nmark(iC+1)-1))+VwormDir;
    for iS=Nmark(iC)+1:Nmark(iC+1)
        if (iS==Nmark(iC)+1)
        idAffect=iS+1; %+-dr affects neighbor on the right only
        elseif (iS==Nmark(iC)+2)
        idAffect=(iS:iS+1); %+-dr affects neighbor on the right and itself
        elseif (iS==Nmark(iC+1))
        idAffect=iS-1; %+-dr affects neighbor on the left only
        elseif (iS==Nmark(iC+1)-1)
        idAffect=(iS-1:iS); %+-dr affects neighbor on the left and itself
        else
        idAffect=(iS-1:iS+1); %+-dr affects neighbors on the left and right
        end
        for iXYZ=1:3
        md.var.coord(iS,iXYZ)=md.var.coord(iS,iXYZ)+dr;
%         [md] = getCosAng(md,'idx',(idAffect));
        [md] = getCosAng(md);
        [VwormAlt(idAffect),~] = potential.Vworm(md.var.cosAng(idAffect),Vin_pm);
        VwormDirAlt=potential.VwormDir(md.var.dir(Nmark(iC)+1:Nmark(iC+1),:),md.pm.Vpm.kappaDir);
        f.int_comp.(Mname)(iS,iXYZ)=f.int_comp.(Mname)(iS,iXYZ)-sum(VwormAlt(idAffect)-Vworm(idAffect))/dr-(VwormDirAlt-VwormDir(iC))/dr;
        md.var.coord(iS,iXYZ)=md.var.coord(iS,iXYZ)-dr;
        VwormAlt(idAffect)=Vworm(idAffect);
        end
    end
end
%----------------------------------------------------------------------------------------
f.int_rand.(Mname)=(rand(md.var.n_coord,3)-0.5)*10000;
% f.int_rand.(Mname)=zeros(md.var.n_coord,3);
%----------------------------------------------------------------------------------------
%% sum up into total force 
f.int_V.(Mname)=Vtot;
f.int_tot.(Mname)=cell(1,1);
f.int_tot.(Mname){1}=f.int_const.(Mname)+f.int_comp.(Mname)+f.int_rand.(Mname);