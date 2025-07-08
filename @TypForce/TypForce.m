classdef TypForce
%==========================================================================
%--------------------------------------------------------------------------
        % TypForce is for quantifying the mechanical interaction between
        % modules like @ModMembrane and @ModSubstrate
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
%==========================================================================
    properties
        int_comp %computed forces
        int_rand %random forces
        int_V %potential
        handle %name handle
        int_const %constant forces, i.e. from linear interpretation potential
        int_stored %stored forces avoid repeating computation
        int_field %field forces
        int_tot %total forces
        n_int %number of forces
        int_avail %availability of a given force
        dt
        otherInfo
        pm %parameters like spring constants
    end
%==========================================================================
%==========================================================================
    methods
%==========================================================================
        function obj = TypForce(M,varargin)
%--------------------------------------------------------------------------
        % TypForce is the initiation function for @TypForce to setup
        % parameters and variables, naming and orgnizing the modules 
        % input: 
        % M - module included such as @ModMembrane
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------            
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('M', @(x) isa(x,'model'));
            ip.parse(M,varargin{:});
%--------------------------------------------------------------------------
            unit=M.unit;
            if isempty(unit)
                unit=nat_unit('erg',ComUnit.nm_to_cm(1),300);
                warning('no unit assigned, using 1nm and 300T as natural units')
            end
%-------------------------------------------------------------------------- 
            [obj] = SetParameter(obj,'unit',unit);
%--------------------------------------------------------------------------
            TypName='TypForce';
%--------------------------------------------------------------------------
            int_path=[M.mod_path filesep '@' TypName];
            D=dir(int_path);
            n_D = numel(D);
            n_int=0;
            is_Mod=false(n_D,1);
            for i=1:n_D
                if numel(D(i).name) > 4
                if strcmp(D(i).name(1:3),'Mod')
                    n_int=n_int+1;
                    is_Mod(i)=true;
                end
                end
            end
            n_int=0;
            for i=1:n_D
                if is_Mod(i)==true
                    n_tem=numel(D(i).name);
                    for i_tem=1:n_tem
                       if strcmp(D(i).name(i_tem),'.') 
                           idx=i_tem;
                           break;
                       elseif strcmp(D(i).name(i_tem),'8')
                           idx=[];
                           break;
                       end
                    end
                    if ~isempty(idx)
                        n_int=n_int+1;
                        obj.int_avail{n_int,1}=D(i).name(1:idx-1);
                    end
                end
            end
            
            if isempty(obj.int_avail{1})
                error(['No force found in @' TypName]);
            else
                obj.int_avail=sort(obj.int_avail);
                obj.int_avail=unique(obj.int_avail);
            end
%--------------------------------------------------------------------------
            field_tem = cell(0,0);
            for i=1:M.n_typ
                for j=1:M.n_int(i)+M.n_mod
                    if strcmp(M.mat_int{i}{j}.type,TypName)
                        field_tem=[field_tem;M.mat_int{i}{j}.name];
                        obj.handle=[obj.handle;{M.mat_int{i}{j}.handle}];
                    end
                end
            end
            obj.n_int=numel(field_tem);
%--------------------------------------------------------------------------
            id_tem=true(obj.n_int,1);
            for i=1:obj.n_int
                int_found=false;
                for j=1:numel(obj.int_avail)
                    if strcmp(field_tem{i},obj.int_avail{j})
                        int_found=true;
                        break;
                    end
                end
                if int_found==false
                    warning([field_tem{i} ' not found in ' TypName]);
                    id_tem(i)=false;
                end
            end
            field_tem=field_tem(id_tem);
            obj.n_int=numel(field_tem);
            obj.handle=obj.handle(id_tem);
%--------------------------------------------------------------------------           
            obj.int_comp=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_rand=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_V=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_const=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_stored=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_field=cell(M.n_mod,1);
            obj.int_tot=cell2struct(cell(obj.n_int,1), field_tem);
            obj.otherInfo=cell2struct(cell(obj.n_int,1), field_tem);
%--------------------------------------------------------------------------
        end
%==========================================================================
        function [obj] = comp_int(obj, M,i_typ, varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj', @(x) isa(x,'TypForce'));
            ip.addRequired('M', @(x) isa(x,'model'));
            ip.addRequired('i_typ', @(x) isnumeric(x));
            ip.parse(obj, M,i_typ,varargin{:});
%--------------------------------------------------------------------------
            for i_int = 1:M.TypForce.n_int
                f_tem = M.TypForce.handle{i_int};
                n_mod=numel(M.mat_int{i_typ}{i_int}.i_mod);
                if n_mod==1
                    mod_tem=M.mod(M.mat_int{i_typ}{i_int}.i_mod);
                    obj = f_tem(obj,mod_tem{:},'mex',Mex);
                elseif n_mod==2
                    obj = interaction(obj,M,M.mat_int{i_typ}{i_int}.name);
                else
                    error('>2 or <1 mod for TypForce' );
                end
            end
        end
%==========================================================================
        %[f,loc_relaxed,abnormal_length,m] = ModMembrane(f,m, varargin);
%         [f,V] = ModClathrin(f,mod, varargin);
        [m] = ModMembrane8Helfrich(f,m, varargin);   
        [f] = var_dt(f,M,varargin);
        [f] = ModSubstrate(f,m, varargin);
        [f] = interaction(f,M,int_name,varargin);
        [f,V] = ModClathrin_ModMembrane(f,M,int_name,varargin);
        [f,V_tot,otherInfo] = ModClathrin_ModFreeParticle(f,M,varargin);
        [id_reps_m_c_leg] = ModClathrin_ModMembrane8getIDrepulsive(f,M,int_name,varargin);
        [id_reps_m_leg] = ModClathrin_ModMembrane8getIDrepSingle(f,M,coord_new,O,varargin);
        [m] = ModMembrane8uK(f,m, varargin);
        [f,V] = ModFreeParticle_ModMembrane(f,M,varargin);
        [f] = storedForce(f,Vname,Vpm,Mname,rPara,varargin);
        [f,Vconst] = constForce(f,Mname,md,idPair,idCoord,varargin);
        [f,Vtot] = ModPolymerChain(f,M,varargin);
%==========================================================================
    end
%==========================================================================
%==========================================================================
methods(Static)
    [V] = ModMembrane8potential1D(r,Vpm,Vcase);
    [V] = ModClathrin8V(M,int_name, varargin);
    [V] = ModClathrin_ModMembrane8V(M,int_name,id_td_mem,varargin);
    [V] = ModClathrin_ModMemAdapter_ModMembrane8V(M,int_name,id_td_mem,varargin);
    [id_td_mem] = ModClathrin_ModMembrane8getADP(M,int_name,varargin);
    [S] = ModRingmem8compS(r1,r2,z1,z2);
    %% Mex
    [f_b,f_p,u_K,kH,f_AV,V] = ModMembraneMex(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,dens);
    [f_b,f_p,u_K,kH,f_AV,V] = ModMembraneAltMex(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,dens);
end
end

