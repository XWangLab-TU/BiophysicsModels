classdef TypChemistry
%==========================================================================
%--------------------------------------------------------------------------
        % TypChemistry is under development for incorporating dynamical
        % systems like signaling pathway simulation
        %   See also model, TypForce
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
%==========================================================================
    properties
        int_comp
        handle
        int_const
        int_stored
        int_tot
        n_int
        int_avail
        dt
        pm
    end
%==========================================================================
%==========================================================================
    methods
%==========================================================================
        function obj = TypChemistry(mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.parse(mod,varargin{:});
%--------------------------------------------------------------------------            
            obj.pm.k_ModMembrane_ModSubstrate=1;
            TypName='TypChemistry';
%--------------------------------------------------------------------------
            int_path=[mod.mod_path filesep '@' TypName];
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
            for i=1:mod.n_typ
                for j=1:mod.n_int(i)+mod.n_mod
                    if strcmp(mod.mat_int{i}{j}.type,TypName)
                        field_tem=[field_tem;mod.mat_int{i}{j}.name];
                        obj.handle=[obj.handle;{mod.mat_int{i}{j}.handle}];
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
            obj.int_const=cell2struct(cell(obj.n_int,1), field_tem);
            obj.int_stored=cell2struct(cell(obj.n_int,1), field_tem);
%--------------------------------------------------------------------------
            obj.int_tot=zeros(mod.n_coord,3);
        end
%==========================================================================
        function [obj] = comp_int(obj, mod,i_typ, varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj', @(x) isa(x,'TypChemistry'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addRequired('i_typ', @(x) isnumeric(x));
            ip.parse(obj, mod,i_typ,varargin{:});
%--------------------------------------------------------------------------
            obj.int_tot=obj.int_tot*0;
            for i_int = 1:mod.TypForce.n_int
                f_tem = mod.TypForce.handle{i_int};
                n_mod=numel(mod.mat_int{i_typ}{i_int}.i_mod);
                if n_mod==1
                    mod_tem=mod.mod(mod.mat_int{i_typ}{i_int}.i_mod);
                    obj = f_tem(obj,mod_tem{:},'mex',Mex);
                elseif n_mod==2
                    obj = interaction(obj,mod,mod.mat_int{i_typ}{i_int}.i_mod,mod.mat_int{i_typ}{i_int}.name);
                else
                    error('>2 or <1 mod for TypForce' );
                end
            end
        end
%==========================================================================        
        [ch] = interaction(ch,id_int,mod);   
%--------------------------------------------------------------------------        
        [mod,changed] = ModClathrin(ch,mod,varargin);
        [mod] = ModClathrin8link(ch,mod,varargin);
        [mod,changed] = ModClathrin8init(ch,mod,r_init,varargin);
%==========================================================================
    end
%==========================================================================
%==========================================================================
methods(Static)
    [waist] = ModClathrin8getWaist(mod,varargin);
    [id_ring_ver1,id_ring_edg1,n_ring_ver1,n_ring_edg1]=ModMembrane8RingMex(edge_all,edge_ctr);
end
end

