function [f,m,Vtot] = ModMembrane(f,mod,varargin)
%--------------------------------------------------------------------------
        % ModMembrane performs the computation of the forces involved in
        % @ModMembrane, including bending force, internal force, pressure
        % force and tension force
        % input: 
        % f - @TypForce
        % mod - @model object
        % output:
        % abnormal_length - whether too long or too short edge appear
        % m - @ModMembrane object with updated force data
        % Vtot - total potential value of m
        % optional:
        % see variable arguments
        %   See also ModMembrane_ModSubstrate
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('Pflat', false, @islogical);
ip.addParameter('split', [], @isstruct);
ip.addParameter('merge', [], @isstruct);
ip.parse(f,mod,varargin{:});
%----------------------------------------------------------------------------------------
m=mod.mod{mod.i_mod.ModMembrane};
%----------------------------------------------------------------------------------------
mex_avail=true;
%----------------------------------------------------------------------------------------
Vpm=m.pm.Vdh;
Vname='VinMembrane2hill';
%----------------------------------------------------------------------------------------
% new method (generalized constant force)
Mname='ModMembrane';
if isempty(f.int_stored.ModMembrane)
    %define parameters for double hill potential
    Vin_pm(1)=Vpm.V0;
    Vin_pm(2)=Vpm.r_1;Vin_pm(3)=Vpm.r_2;
    Vin_pm(4)=Vpm.k_w;
    Vin_pm(5)=Vpm.e_b1;Vin_pm(6)=Vpm.e_b2;
    Vin_pm(7)=Vpm.e_w;
    Vin_pm(8)=Vpm.k_b11; Vin_pm(9)=Vpm.k_b12; Vin_pm(10)=Vpm.k_b21;Vin_pm(11)=Vpm.k_b22;
    Vin_pm(12)=Vpm.rb_11;Vin_pm(13)=Vpm.rb_12;Vin_pm(14)=Vpm.rb_21;Vin_pm(15)=Vpm.rb_22;
    f=storedForce(f,Vname,Vin_pm,Mname,[Vpm.r_1 m.pm.dr Vpm.r_2],'rsq_std', m.pm.f_const_rsq_std,'std_std', m.pm.f_const_std_std);
end
idPair=m.var.edge_all(m.var.id_on_edg,:);
idCoord=m.var.id_on_coord;
[f,Vtot] = constForce(f,Mname,m,idPair,idCoord);
abnormal_length=nan; %no need for abnormal_length

%old method
% if isempty(f.int_stored.ModMembrane)
% %     f.int_stored.ModMembrane = ModMembrane8Const(f,m);
% end
% d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
% r = sqrt(sum(d.^2,2));
% u = d./r;
% i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
% i = floor(r/m.pm.dr+0.5)-i_shift;
% 
% f_edg=f.int_stored.ModMembrane.fn(i).*u;
% f.int_const.ModMembrane=zeros(m.var.n_coord,3);
% Vtot=sum(f.int_stored.ModMembrane.Vn(i));
% for i_coord = 1:numel(m.var.id_on_coord)
%     f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)-sum(f_edg(m.var.edge_all(m.var.id_on_edg,1)==m.var.id_on_coord(i_coord),:),1);
%     f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)+sum(f_edg(m.var.edge_all(m.var.id_on_edg,2)==m.var.id_on_coord(i_coord),:),1);
% end
%----------------------------------------------------------------------------------------
% id_tem1 = r<Vpm.rl_min;  id_tem2 = r>Vpm.rl_max; 
% n_tem1 = length(id_tem1(id_tem1));
% n_tem2 = length(id_tem2(id_tem2));
% if (n_tem1 == 0) && (n_tem2 == 0)
%     abnormal_length=false;
% else
%     abnormal_length=true;
% end
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
% f.int_rand.ModMembrane = randn(m.var.n_coord,3)*sqrt(2*m.pm.mu*m.pm.kBT*m.pm.dt)/m.pm.dt/m.pm.mu; %XXX
f.int_rand.ModMembrane = zeros(m.var.n_coord,3);

if mex_avail == false
%%
[m] = f.ModMembrane8Helfrich(m,'idx',m.var.id_on_coord);
coord_org = m.var.coord;
H_org = m.var.f.H;
Vtot=Vtot+H_org;
K_org=m.var.f.K;
kH_org = m.var.f.kH;
A_org = m.var.f.A;
f.int_comp.ModMembrane = zeros(m.var.n_coord,3);
m.var.f.dH = m.var.f.dH*0;
dr_H=0.00001;
for i_on = 1:m.var.n_on_coord
    i=m.var.id_on_coord(i_on);
    for k=1:3
       m.var.coord = coord_org; m.var.coord(i,k) = m.var.coord(i,k) + dr_H;
       idx = [i,m.var.j_T(i,1:m.var.val(i))];
       [~,id_tem]=intersect(idx,m.var.id_on_coord);
       idx=idx(id_tem);
       m = f.ModMembrane8Helfrich(m, 'idx',idx,'init',false);
       H_new = m.var.f.H;
%        id_tem=1:m.var.n_coord;
%        for i_tem=1:max(size(idx))
%            id_tem(id_tem==idx(i_tem))=[];
%        end
%        norm(kH_org(id_tem,:) - m.var.f.kH(id_tem,:))
%        pause
       f.int_comp.ModMembrane(i,k) = -m.pm.k_c*(H_new-H_org)/dr_H;
       m.var.f.dH(i,k) = (H_new-H_org)/dr_H;
       m.var.f.kH(idx)=kH_org(idx);
       m.var.f.A(idx)=A_org(idx);
       m.var.f.K(idx,:)=K_org(idx,:);
    end
end
%%
K_tem = m.var.f.K./m.var.f.A;
K_n=sqrt(sum(K_tem.*K_tem,2));
u_K = K_tem./K_n;

m.var.coord = coord_org;
else
    pmc=zeros(10,1);
    pmc(1) = m.pm.nt;
    pmc(2) = m.pm.dt;
    pmc(3) = m.pm.P;
    pmc(4) = m.pm.k_c;
    pmc(5) = m.pm.k_e;
    pmc(6) = m.pm.dr;
    pmc(7) = m.pm.k_V;
    pmc(8) = m.pm.k_A;
    pmc(9) = m.pm.V0;
    pmc(10) = m.pm.A0;
    pmc(11) = m.pm.nAVmean;
    pmc(12) = m.pm.k_a;
       
    j_T = m.var.j_T; j_T(isnan(j_T)) = 0;
    
%     V_A_Vol: V, membrane free energy, A, Area, Vol, volume
    if m.pm.GaussianCurv==false
       [f_b,f_p,u_K,kH,f_AV,V_A_Vol]=f.ModMembraneMex...
       (m.var.coord,...
        pmc,...
        m.var.edge_all,...
        m.var.face_unq,...
        j_T,...
        m.var.id_on_coord,...
        m.var.n_node',...
        m.var.T_s',...
        m.var.T_e',...
        m.var.dens);
    else
       [f_b,f_p,u_K,kH,f_AV,V_A_Vol]=f.ModMembraneAltMex...
       (m.var.coord,...
        pmc,...
        m.var.edge_all,...
        m.var.face_unq,...
        j_T,...
        m.var.id_on_coord,...
        m.var.n_node',...
        m.var.T_s',...
        m.var.T_e',...
        m.var.dens);
    end
    %%
%     fig=figure('units','normalized','outerposition',[0 0 1 1]);
%     subplot(1,2,1);
%     plot(m,'f',fig,'LineStyle','-','col',m.var.dens/max(m.var.dens));
%     subplot(1,2,2);
%     plot(m,'f',fig,'LineStyle','-','col',sum(f_AV.^2,2)/max(sum(f_AV.^2,2)));   
%     quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),f_AV(:,1),f_AV(:,2),f_AV(:,3)); hold on;
%%
    m.var.f.kH=kH;
    m.var.f.u_K=u_K;
    m.var.f.f_AV=f_AV;
    if ip.Results.Pflat==true
    f_p(m.var.coord(:,3)<0.0001,:)=0;
    end
    f.int_comp.ModMembrane=f_b+f_p+f_AV;
    Vtot=Vtot+V_A_Vol(1); 
    m.var.f.A=V_A_Vol(2);
    m.var.f.V=V_A_Vol(3);
    %pure bending energy: m.pm.k_c*sum(0.5*(2*m.var.f.kH).^2.*m.var.f.A)
end
% extra forces to prevent unsplitable and unmergeable edges
%%
% try
% if m.pm.doRemeshExtraCtrl==true
%     m.pm.remeshExtraCtrl=5;
%     idEdgAll=1:m.var.n_edg;
%     CanSplitMerge=logical(m.var.CanSplitMerge);
%     idEdgNoSplit=idEdgAll(~CanSplitMerge(:,1))';
%     idEdgNoMerge=idEdgAll(~CanSplitMerge(:,2))';
%     nNoSplit=numel(idEdgNoSplit);
%     nNoMerge=numel(idEdgNoMerge);
% %     disp(min(m.var.CanSplitMerge));
%     if nNoMerge>0
%         r=(m.var.coord(m.var.edge_all(idEdgNoMerge,2),:)-m.var.coord(m.var.edge_all(idEdgNoMerge,1),:));
%         r=sqrt(sum(r.^2))
%         idNoMerge=(1:nNoMerge)';
%         idNoMerge(r>m.pm.l0)=[];
%         idEdgNoMerge=idEdgNoMerge(idNoMerge);
%         idVerNoMerge=unique(m.var.edge_all(idEdgNoMerge,:));
%         f.int_comp.ModMembrane(idVerNoMerge,:)=f.int_comp.ModMembrane(idVerNoMerge,:)+...
%                                                m.pm.remeshExtraCtrl*f.int_const.ModMembrane(idVerNoMerge,:);
%     end
%     if nNoSplit>0
%         r=(m.var.coord(m.var.edge_all(idEdgNoSplit,2),:)-m.var.coord(m.var.edge_all(idEdgNoSplit,1),:));
%         r=sqrt(sum(r.^2));
%         idNoSplit=(1:nNoSplit)';
%         idNoSplit(r<m.pm.l0)=[];
%         idEdgNoSplit=idEdgNoSplit(idNoSplit);
%         idVerNoSplit=unique(m.var.edge_all(idEdgNoSplit,:));
%         f.int_comp.ModMembrane(idVerNoSplit,:)=f.int_comp.ModMembrane(idVerNoSplit,:)+...
%                                                m.pm.remeshExtraCtrl*f.int_const.ModMembrane(idVerNoSplit,:);
%     end
% end
% catch
%     disp('wrong');
%     disp('wrong');
% end
%%
f.int_V.ModMembrane=Vtot;
f.int_tot.ModMembrane=cell(1,1);
%for converting to CS-based remeshing
f.int_tot.ModMembrane{1}=f.int_const.ModMembrane*abs(m.pm.remeshScheme-1)+f.int_comp.ModMembrane+f.int_rand.ModMembrane; 