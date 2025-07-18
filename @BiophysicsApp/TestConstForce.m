function [] = TestConstForce(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', 20000, @isnumeric); %simulation steps for membrane relaxation
ip.addParameter('nRep', 12, @isnumeric); %repeat number
ip.addParameter('xyzLim', [-15 15;-15 15;-15 15], @isnumeric); %figure axis range
ip.addParameter('viewAng', [-50 10], @isnumeric); %figure view angle
ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = 'D:\Matlab_codes\data\computation';
% dirMod='\\lamella.biohpc.swmed.edu\s171152\codes\matlab\mine\git\DanuserLab\biophysicsmodels';addpath(dirMod);
%--------------------------------------------------------------------------
%%
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); %baseline length 1000;
m=ModMembrane(true,3,0,'unit',u);
pc=ModPolymerChain('dyn2','unit',u,'nChain',1,'nSubunit',20); 
p=ModFreeParticle(1,'unit',u,'coordAssigned',pc.var.coord(1,:));
b=ModBrownianMotor('flash',1,'ModFreeParticle','ModPolymerChain','dyn2','unit',u);
%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{p,b},{pc,b}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],m,pc,p,b);
%--------------------------------------------------------------------------
%%
f=M.TypForce;
Vname='VinMembrane2hill';
Mname='ModMembrane';
pm=M.mod{M.i_mod.ModMembrane}.pm;
Vpm(1)=pm.Vdh.V0;
Vpm(2)=pm.Vdh.r_1;Vpm(3)=pm.Vdh.r_2;
Vpm(4)=pm.Vdh.k_w;
Vpm(5)=pm.Vdh.e_b1;Vpm(6)=pm.Vdh.e_b2;
Vpm(7)=pm.Vdh.e_w;
Vpm(8)=pm.Vdh.k_b11; Vpm(9)=pm.Vdh.k_b12; Vpm(10)=pm.Vdh.k_b21;Vpm(11)=pm.Vdh.k_b22;
Vpm(12)=pm.Vdh.rb_11;Vpm(13)=pm.Vdh.rb_12;Vpm(14)=pm.Vdh.rb_21;Vpm(15)=pm.Vdh.rb_22;
f=storedForce(f,Vname,Vpm,Mname,[pm.Vdh.r_1 pm.dr pm.Vdh.r_2],'plot_or_not',true,'rsq_std', pm.f_const_rsq_std,'std_std', pm.f_const_std_std);
m=M.mod{M.i_mod.ModMembrane};
[fc] = ModMembrane8Const(f,m, 'plot_or_not',true);
%--------------------------------------------------------------------------
f_save=f;
f.int_stored.ModMembrane = fc;
d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
u = d./r;
i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
i = floor(r/m.pm.dr+0.5)-i_shift;

f_edg=f.int_stored.ModMembrane.fn(i).*u;
f.int_const.ModMembrane=zeros(m.var.n_coord,3);
Vtot=sum(f.int_stored.ModMembrane.Vn(i));
for i_coord = 1:numel(m.var.id_on_coord)
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)-sum(f_edg(m.var.edge_all(m.var.id_on_edg,1)==m.var.id_on_coord(i_coord),:),1);
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)+sum(f_edg(m.var.edge_all(m.var.id_on_edg,2)==m.var.id_on_coord(i_coord),:),1);
end
f_alt=f;

f=f_save;
idPair=m.var.edge_all(m.var.id_on_edg,:);
idCoord=m.var.id_on_coord;
md=m;
[f,Vconst] = constForce(f,Mname,md,idPair,idCoord);
norm(f.int_const.ModMembrane-f_alt.int_const.ModMembrane)
%----------------------------------------------------------------------------------------
%%

%----------------------------------------------------------------------------------------
end