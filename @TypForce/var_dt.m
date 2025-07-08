function [f] = var_dt(f,m,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addParameter('f_const_only', false, @islogical);
ip.addParameter('case_f', 1, @isnumeric);
ip.addParameter('iMethod', 1, @isnumeric);
ip.parse(f,m,varargin{:});
%--------------------------------------------------------------------------
f_const_only=ip.Results.f_const_only;
case_f=ip.Results.case_f;
iMethod=ip.Results.iMethod;
%--------------------------------------------------------------------------
%========================================================================================================================== compute dt
%%
%--------------------------------------------------------------------------------------------------------membrane
    d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
    r = sqrt(sum(d.^2,2));
    i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
    i = floor(r/m.pm.dr+0.5)-i_shift;
    r=floor(r/m.pm.dr+0.5)*m.pm.dr;
    if f_const_only==true
        f_mod=f.int_const.ModMembrane;
    else
        if case_f==1
            f_mod=f.int_const.ModMembrane+f.int_comp.ModMembrane+f.int_rand.ModMembrane+f.int_comp.ModMembrane_ModSubstrate{1};
        elseif case_f==2
            f_mod=f.int_const.ModMembrane+f.int_comp.ModMembrane+f.int_rand.ModMembrane...
                 +f.int_comp.ModClathrin_ModMembrane{2};
        elseif case_f==3
            f_mod=f.int_const.ModMembrane+f.int_comp.ModMembrane+f.int_rand.ModMembrane...
                 +f.int_comp.ModClathrin_ModMemAdapter_ModMembrane{2};
        elseif case_f==4
            f_mod=f.int_const.ModMembrane+f.int_rand.ModMembrane+f.int_comp.ModMemAdapter_ModMembrane;
        elseif case_f==5
             f_mod=f.int_const.ModMembrane+f.int_rand.ModMemAdapter_ModMembrane+f.int_comp.ModMemAdapter_ModMembrane...
                 +f.int_comp.ModClathrin_ModMemAdapter_ModMembrane{2};
        elseif case_f==6
             f_mod=f.int_const.ModMembrane+f.int_rand.ModMembrane+f.int_comp.ModMemAdapter_ModMembrane...
                 +f.int_comp.ModClathrin_ModMemAdapter_ModMembrane{2}+f.int_comp.ModMembrane_ModSubstrate{1};     
        end
    end
    f1=f_mod(m.var.edge_all(m.var.id_on_edg,1),:);
    f2=f_mod(m.var.edge_all(m.var.id_on_edg,2),:);
    df=(f2-f1);
    SorE=sum(df.*d,2);  %Shrink or Extend
    iSorE=true(size(SorE)); %true: Shrink false: Extend
    iSorE(SorE>0)=false;
    iSorE(SorE<0)=true;
    b=sum(2*m.pm.mu*d.*df,2);
    if iMethod==1
        c_1=r.^2-(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-2)+f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-1))).^2;
        c_2=r.^2-(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i))  +f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)+1))).^2;
        %c_1:shrink, c_2: extend
        c=c_1;
        c(iSorE==false)=c_2(iSorE==false);
        a=m.pm.mu^2*sum(df.^2,2);
        Delta=b.^2-4*a.*c;
        idWrong=Delta<0;
        c(idWrong & (iSorE==false))=c_1(idWrong & (iSorE==false));
        c(idWrong & (iSorE==true))=c_2(idWrong & (iSorE==true));
        Delta(idWrong)=b(idWrong).^2-4*a(idWrong).*c(idWrong);
        dt=[(-b+sqrt(Delta))./(2*a),(-b-sqrt(Delta))./(2*a)];
        dt(dt<0)=inf;
        dt=min(dt,[],2);
        dt=min(dt);
        %% plot check for solution of Delta
%         figure;
%         [~,in]=min(Delta);
%         x2=m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:);hold on;
%         x1=m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:);hold on;
%         x2t=x2(in,:);
%         x1t=x1(in,:);
%         dtTem=0.00001;  
%         scatter3(x2(in,1),x2(in,2),x2(in,3),'filled','markerfacecolor',[1 0 0]);hold on;
%         scatter3(x1(in,1),x1(in,2),x1(in,3),'filled','markerfacecolor',[0 0 1]);hold on;
%         r2Tem=[];
%         for it=1:1000
%            x2t=x2t+dtTem*m.pm.mu*f2(in,:);
%            x1t=x1t+dtTem*m.pm.mu*f1(in,:);
%            scatter3(x2t(1),x2t(2),x2t(3),'filled','markerfacecolor',[1 0 0]);hold on;
%            scatter3(x1t(1),x1t(2),x1t(3),'filled','markerfacecolor',[0 0 1]);hold on;
%            r2Tem=[r2Tem;r(in).^2-(norm(x2t-x1t))^2];
%         end
%         quiver3(x2(in,1),x2(in,2),x2(in,3),f2(in,1),f2(in,2),f2(in,3));hold on;
%         quiver3(x1(in,1),x1(in,2),x1(in,3),f1(in,1),f1(in,2),f1(in,3));hold on;
%         figure;
%         plot(r2Tem,'.');
        %%        
    else
        %%
        c_1=(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-2)+f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)-1))).^2-r.^2;
        c_2=(0.5*(f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i))  +f.int_stored.ModMembrane.rg(f.int_stored.ModMembrane.in(i)+1))).^2-r.^2;
        dt_det=[c_1 c_2]./b;
        id_keep=dt_det>0;
        dt_final=sum(dt_det.*id_keep,2);    
        dt=min(dt_final);
    end
%--------------------------------------------------------------------------------------------------------
f.dt=dt;