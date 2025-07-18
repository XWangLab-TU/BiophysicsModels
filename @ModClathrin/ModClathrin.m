classdef ModClathrin
%==========================================================================
%==========================================================================     
    properties
       var
       pm
       prop
       failInfo
    end
%==========================================================================
%==========================================================================     
    methods
        function obj = ModClathrin(varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addParameter('unit', [], @isobject);
            ip.parse(varargin{:});
%--------------------------------------------------------------------------            
            unit=ip.Results.unit;
            if isempty(unit)
                unit=ComUnit('erg',ComUnit.nm_to_cm(1),300,ComUnit.kBT_to_erg(1,300)); 
                warning('no unit assigned, using 1nm, 1kBT and 300K as natural units');
            end
            kBT=unit.unit.kBT; %(in erg,K)
%--------------------------------------------------------------------------
            obj.prop = {'RigidBody','needBreakOff'};
            obj.failInfo='';
            obj.pm=struct('mu_r',[],'mu_t',[],'eta',0.01,'k_on',10000,... %eta:Solvent viscosity (poise or g⋅cm−1⋅s−1)
                          'O_test_all',[],'n_test_all',[],'ang_a_test_all',[],'a_test_all',[],...
                          'r_min',20,...
                          'cos_min',0.9,...
                          'l0',1.5,...
                          'l1',2.5,...
                          'r_min_for_link',25,... %25nm
                          'n_init_adp',3,...
                          'n_min_adp',3,...
                          'n_max_adp',3,...
                          'psi',-0.4... %11 degree(default):-0.191986217719376 
                          ); 
            % mu_r=Dr/kBT,[Dr]=s^-1 
%             Dr=[1.647E+04 -8.673E-01 -5.987E-02;...
%                 -8.666E-01  1.647E+04  5.337E-02;...
%                 -6.023E-02  5.280E-02  1.024E+04]; %default

            Dr=[2.288E+04 -1.052E-01  2.055E-01;...
                -1.083E-01  2.288E+04  9.004E-01;...
                2.049E-01  8.983E-01  1.605E+04 ];   %psi=-0.4
         
            obj.pm.mu_r=Dr/kBT;
            % mu_t=Dt/kBT,[Dr]=s^-1
%             Dt=[1.337E-07 -3.859E-12  1.577E-13;...
%                 -3.845E-12  1.337E-07 -1.698E-13;...
%                 1.588E-13 -1.690E-13  1.091E-07]; %default
           Dt=[1.296E-07 -3.586E-14  4.578E-13;...     
               -6.275E-14  1.296E-07  5.190E-13;...    
                4.620E-13  5.247E-13  1.313E-07];  %psi=-0.4
            obj.pm.mu_t=Dt/kBT;
           %%
            d_ang=0.1;
            dPhi=d_ang*pi;
            dtheta=d_ang*pi;
            dphi=d_ang*pi;
            P_grid=0:dPhi:pi;
            theta_grid=0:dtheta:1*pi;
            phi_grid=0:dphi:2*pi;
            obj.pm.ang_a_test_all=ComMath.grid3(P_grid',theta_grid',phi_grid');
            obj.pm.n_test_all=size(obj.pm.ang_a_test_all,1);
            obj.pm.a_test_all=obj.get_a_from_ang(obj.pm.ang_a_test_all(:,1),obj.pm.ang_a_test_all(:,2),obj.pm.ang_a_test_all(:,3));
            obj.pm.O_test_all=obj.Omega(obj.pm.a_test_all,obj.pm.ang_a_test_all(:,1));
            u=zeros(obj.pm.n_test_all,3);
            u(:,3)=1;
            for i=1:obj.pm.n_test_all
                u(i,:)=obj.get_r_from_a(u(i,:),1,obj.pm.O_test_all(:,:,i));
            end 
            d_xyz=0.2;
            range=cell(3,1);
            range{1}=(-1:d_xyz:1)';range{2}=(-1:d_xyz:1)';range{3}=(-1:d_xyz:1)';
            coord=ComMath.grid3(range{1},range{2},range{3});
            [id_u] = ComMath.getMeshID(u,coord,range,d_xyz);
            id_not_nan=~isnan(id_u);
            u=u(id_not_nan,:);
            id_u=id_u(id_not_nan);
            [~,ia,~] = unique(id_u);
            %u_unq = u(ia,:);scatter3(u_unq(:,1),u_unq(:,2),u_unq(:,3),'filled');
            
            obj.pm.ang_a_test_all=obj.pm.ang_a_test_all(id_not_nan,:);
            obj.pm.ang_a_test_all=obj.pm.ang_a_test_all(ia,:);
            obj.pm.n_test_all=size(obj.pm.ang_a_test_all,1);
            obj.pm.a_test_all=obj.get_a_from_ang(obj.pm.ang_a_test_all(:,1),obj.pm.ang_a_test_all(:,2),obj.pm.ang_a_test_all(:,3));
            obj.pm.O_test_all=obj.Omega(obj.pm.a_test_all,obj.pm.ang_a_test_all(:,1));
%             u=zeros(obj.pm.n_test_all,3);
%             u(:,3)=1;
%             for i=1:obj.pm.n_test_all
%                 u(i,:)=obj.get_r_from_a(u(i,:),1,obj.pm.O_test_all(:,:,i));
%             end 
%             scatter3(u(:,1),u(:,2),u(:,3),'filled');
%==========================================================================            
            %%
            l0=obj.pm.l0; l1=obj.pm.l1; 
            %phi=-2*pi/3; psi=0/90*0.5*pi;  Xi= -3*(0.5*pi+psi - pi/6)/2*0.;
            phi=-2*pi/3; psi=obj.pm.psi;  Xi= -3*(0.5*pi+psi - pi/6)/2*1;
            coord=cell(10);
            coord_std=[(0:8)'*2*l0,zeros(9,1),zeros(9,1)];
            phi0=2*pi/3;
            coord{1}=coord_std;
            coord{2}=obj.rotation(coord{1},0,phi0);
            coord{3}=obj.rotation(coord{2},0,phi0);
            coord{4}=-coord{1}; 
            coord{4}=obj.rotation(coord{4},0,phi);
            coord{4}=coord{4}+coord{1}(end,:); 
            coord{5}=-coord{2}; 
            coord{5}=obj.rotation(coord{5},0,phi);
            coord{5}=coord{5}+coord{2}(end,:); 
            coord{6}=-coord{3}; 
            coord{6}=obj.rotation(coord{6},0,phi);
            coord{6}=coord{6}+coord{3}(end,:); 
            phi1=-2*pi/3;
            coord{7}=-(coord{4}(end,:)-coord{4}(end-1,:))/sqrt(sum((coord{4}(end,:)-coord{4}(end-1,:)).^2,2))*(l0+l1); 
            coord{7}=obj.rotation(coord{7},0,phi1);
            coord{7}=coord{7}+coord{4}(end,:);
            coord{8}=-(coord{5}(end,:)-coord{5}(end-1,:))/sqrt(sum((coord{5}(end,:)-coord{5}(end-1,:)).^2,2))*(l0+l1); 
            coord{8}=obj.rotation(coord{8},0,phi1);
            coord{8}=coord{8}+coord{5}(end,:);
            coord{9}=-(coord{6}(end,:)-coord{6}(end-1,:))/sqrt(sum((coord{6}(end,:)-coord{6}(end-1,:)).^2,2))*(l0+l1); 
            coord{9}=obj.rotation(coord{9},0,phi1);
            coord{9}=coord{9}+coord{6}(end,:);
            coord{10}=[0 0 0];
            for i=1:6
                coord{i}(1,:)=[];
            end

            for i=1:3:7
                coord{i}=obj.rotation(coord{i},0,0.5*pi);
                coord{i}=obj.rotation(coord{i},0.5*pi,0);
            end
            coord{4}=obj.rotation(coord{4},0,Xi);
            coord{7}=obj.rotation(coord{7},0,Xi);
            for i=1:3:7
                coord{i}=obj.rotation(coord{i},-0.5*pi,0);
                coord{i}=obj.rotation(coord{i},0,-0.5*pi);
            end
            for i=2:3:8
                coord{i}=obj.rotation(coord{i},0,-2*pi/3+0.5*pi);
                coord{i}=obj.rotation(coord{i},0.5*pi,0);
            end
            coord{5}=obj.rotation(coord{5},0,Xi);
            coord{8}=obj.rotation(coord{8},0,Xi);
            for i=2:3:8
                coord{i}=obj.rotation(coord{i},-0.5*pi,0);
                coord{i}=obj.rotation(coord{i},0,2*pi/3-0.5*pi);
            end
            for i=3:3:9
                coord{i}=obj.rotation(coord{i},0,-4*pi/3+0.5*pi);
                coord{i}=obj.rotation(coord{i},0.5*pi,0);
            end
            coord{6}=obj.rotation(coord{6},0,Xi);
            coord{9}=obj.rotation(coord{9},0,Xi);
            for i=3:3:9
                coord{i}=obj.rotation(coord{i},-0.5*pi,0);
                coord{i}=obj.rotation(coord{i},0,4*pi/3-0.5*pi);
            end
                       
            for i=1:3:7
                coord{i}=obj.rotation(coord{i},0,0.5*pi);
                coord{i}=obj.rotation(coord{i},psi,0);
                coord{i}=obj.rotation(coord{i},0,-0.5*pi);
            end
            for i=2:3:8
                coord{i}=obj.rotation(coord{i},0,-2*pi/3+0.5*pi);
                coord{i}=obj.rotation(coord{i},psi,0);
                coord{i}=obj.rotation(coord{i},0,2*pi/3-0.5*pi);
            end
            for i=3:3:9
                coord{i}=obj.rotation(coord{i},0,-4*pi/3+0.5*pi);
                coord{i}=obj.rotation(coord{i},psi,0);
                coord{i}=obj.rotation(coord{i},0,4*pi/3-0.5*pi);
            end
            
            
            for i=1:10
                if (i>=7)&&(i<=9)
                    coord{i}=[coord{i},l1*ones(size(coord{i},1),1)];
                else
                    coord{i}=[coord{i},l0*ones(size(coord{i},1),1)];
                end
            end
            coord_all=[];
            for i=1:10
                coord_all=[coord_all;coord{i}];
            end
%==========================================================================
            u=unit;
            u=convert(u,'length',ComUnit.nm_to_cm(coord_all));
            coord_all=u.unit_nat.length;
            u=convert_high_order(u,obj.pm.mu_r,0,-1);
            obj.pm.mu_r=u.unit_nat_any;
            u=convert_high_order(u,obj.pm.mu_t,2,-1);
            obj.pm.mu_t=u.unit_nat_any;
            u=convert(u,'length',ComUnit.nm_to_cm(obj.pm.l0));
            obj.pm.l0=u.unit_nat.length;
            u=convert(u,'length',ComUnit.nm_to_cm(obj.pm.l1));
            obj.pm.l1=u.unit_nat.length;
            u=convert(u,'length',ComUnit.nm_to_cm(obj.pm.r_min_for_link));
            obj.pm.r_min_for_link=u.unit_nat.length;
%==========================================================================
           %%
            [obj.var]=setVar(obj,false,'coord_org',coord_all);
%==========================================================================   
        end
%==========================================================================         
        [f] = plot(obj,varargin);
        [var] = setVar(obj,update,varargin);
        [waist] = getWaist(obj,varargin);
        [CtrlPt] = getCtrlPt(obj,varargin);
        [foot] = getFoot(obj,varargin);
        [obj,changed] = link(obj,varargin);
        [T] = getTightness(obj,M,varargin);
%==========================================================================
    end
%==========================================================================    
    methods (Static)
        function mod_name = identify(varargin)
            mod_name = 'ModClathrin';
        end
        [coord_final]=rotation(coord,theta,phi);
        [a]=get_a_from_ang(Phi,theta_a,phi_a);
        [O]=Omega(a,Phi);
        [r]=get_r_from_a(r,n,O);
        [E]=getE(a,Phi);
    end
%==========================================================================    
end

