classdef model
%==========================================================================
%--------------------------------------------------------------------------
        % model stores the data structure of a serious of modules, e.g.
        % @ModMembrane and @ ModSubstrate as cell array, also included
        % interactions of @TypForce and @TypChemistry
        %   See also ModMembrane, ModSubstrate
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------  
%==========================================================================    
    properties
        mod %inividual module
        name 
        changed
        pm
        mod_path
        mod_avail
        int_avail
        Mex
        n_mod %number module
        i_mod %index of module
        n_typ
        i_typ
        n_int
        i_int
        n_elm
        mat_int   %1st idx: type; 2nd idx: interaction; 
        n_coord
        i_coord
        TypForce %@TypForce interaction 
        TypChemistry %@TypChemistry interaction
        unit %@ComUnit as natural unit
        Mesh %cubic mesh in 3D
        t %time
    end
%==========================================================================
%==========================================================================    
    methods
        function obj = model(int_info,dir_mod,unit,lim_xyz,varargin)
%--------------------------------------------------------------------------
        % model is the initiation function for @model to setup
        % parameters and variables, naming and orgnizing the modules 
        % input: 
        % int_info - variable for naming purpose
        % dir_mod - directory of modules included
        % unit - @ComUnit object for defining natural unit
        % lim_xyz - spatial  dimension of cubic mesh
        % varargin:
        % all modules included
        %   See also TypForce, TypChemistry
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('int_info', @(x) isstruct(x));
            ip.addRequired('dir_mod', @(x) ischar(x));
            ip.addRequired('unit', @(x) isa(x,'ComUnit'));
            ip.addRequired('lim_xyz', @(x) isnumeric(x));
            ip.parse(int_info,dir_mod,unit,lim_xyz);
%--------------------------------------------------------------------------
            obj.mod_path = dir_mod;
            obj.Mex=Mex;
            D=dir(dir_mod);	n_D = numel(D);
            obj.t=0;
%--------------------------------------------------------------------------  
           obj.unit=unit;
%-------------------------------------------------------------------------- 
           obj.pm=struct('kBT',1);
%--------------------------------------------------------------------------
           lim_xyz=ip.Results.lim_xyz;
           dr=1;
           lim_x=lim_xyz(1,:);lim_y=lim_xyz(2,:);lim_z=lim_xyz(3,:);          
           obj.Mesh=struct('d',dr,'coord',[],'range',[]);
           rang_x=(lim_x(1):dr:lim_x(2));rang_y=(lim_y(1):dr:lim_y(2));rang_z=(lim_z(1):dr:lim_z(2));
           rang_x=rang_x';rang_y=rang_y';rang_z=rang_z';
           obj.Mesh.coord=ComMath.grid3(rang_x,rang_y,rang_z);
           obj.Mesh.range=cell(3,1);
           obj.Mesh.range{1}=rang_x;
           obj.Mesh.range{2}=rang_y;
           obj.Mesh.range{3}=rang_z;
%--------------------------------------------------------------------------
            n_int=0;
            is_Int=false(n_D,1);
            for i=1:n_D
                if numel(D(i).name) > 4
                if strcmp(D(i).name(2:4),'Typ')
                    n_int=n_int+1; is_Int(i)=true;
                end
                end
            end
            if n_int==0; error('No interaction is found in the given dir_mod'); end
            obj.int_avail = cell(n_int,1);
            n_int=0;
            for i=1:n_D
                if is_Int(i)==true
                    n_int=n_int+1; obj.int_avail{n_int}=D(i).name(2:end);
                end
            end
            obj.int_avail=sort(obj.int_avail);
%--------------------------------------------------------------------------            
            n_mod=0;
            is_Mod=false(n_D,1);
            for i=1:n_D
                if numel(D(i).name) > 4
                if strcmp(D(i).name(2:4),'Mod')
                    n_mod=n_mod+1; is_Mod(i)=true;
                end
                end
            end
            if n_mod==0; error('No module is found in the given dir_mod'); end
            obj.mod_avail = cell(n_mod,1);
            n_mod=0;
            for i=1:n_D
                if is_Mod(i)==true
                    n_mod=n_mod+1; obj.mod_avail{n_mod}=D(i).name(2:end);
                end
            end
            obj.mod_avail=sort(obj.mod_avail);
%--------------------------------------------------------------------------            
            obj.n_mod=numel(varargin);
%--------------------------------------------------------------------------
            obj.n_typ=numel(int_info);
            obj.i_typ=1:obj.n_typ;
%--------------------------------------------------------------------------
            n_coorp_max=0;
            for i=1:obj.n_typ
                for j=1:size(int_info(i).IntList,1)
                   i_tem=size(int_info(i).IntList{j},2);
                   if i_tem>n_coorp_max
                       n_coorp_max=i_tem;
                   end
                end
            end
            
            if n_coorp_max==0
                warning('no coorporative interaction assigned');
                obj.n_elm=[]; 
            else
                obj.n_elm=zeros(obj.n_typ,n_coorp_max);
            end
            obj.n_int=zeros(obj.n_typ,1);
            
            for i=1:obj.n_typ
                obj.n_int(i)=numel(int_info(i).IntList);
                if isempty(int_info(i).IntList{1})
                   obj.n_int(i)=0;
                end
                for j=1:obj.n_int(i)
                    obj.n_elm(i,j)=size(int_info(i).IntList{j},2);
                end
            end  
%--------------------------------------------------------------------------
            obj.changed=false(obj.n_mod,1);
%--------------------------------------------------------------------------            
            mod_list=cell(obj.n_mod,1);
            obj.mod = cell(obj.n_mod,1);
            obj.name=cell(obj.n_mod,1);
            for i=1:obj.n_mod
                mod_list{i} = class(varargin{i});
            end
            [~,idx]=sort(mod_list);
            fprintf('mod_list is sorted alphabetically\n');
            for i=1:obj.n_mod
                obj.mod{i} = varargin{idx(i)};
                obj.name{i}=class(obj.mod{i});
            end          
%--------------------------------------------------------------------------            
            field_tem = cell(0,0);
            for i=1:obj.n_mod
                field_tem=[field_tem;class(obj.mod{i})];
            end
            obj.i_mod=cell2struct(cell(obj.n_mod,1), field_tem);
            for i=1:obj.n_mod
                obj.i_mod.(class(obj.mod{i}))=i;
            end
%--------------------------------------------------------------------------
            obj.mat_int=cell(obj.n_typ,1);
            for i=1:obj.n_typ
                obj.mat_int{i}=cell(obj.n_int(i)+obj.n_mod,1);
            end
%--------------------------------------------------------------------------
            for i=1:obj.n_typ
                for j=1:obj.n_int(i)
                    name_tem=cell(obj.n_elm(i,j),1);
                    for k=1:obj.n_elm(i,j)
                        name_tem{k}=int_info(i).IntList{j}{k};
                    end
                    [name_all,id_mod] = model.getModCoord(name_tem,obj.mod);
                    obj.mat_int{i}{j}=struct('type',int_info(i).TypName,'name',name_all,'handle',str2func(name_all),'i_mod',id_mod);
                end
                for j=obj.n_int(i)+1:obj.n_int(i)+obj.n_mod
                    obj.mat_int{i}{j}=struct('type',int_info(i).TypName,'name',class(obj.mod{j-obj.n_int(i)}),...
                                             'handle',str2func(class(obj.mod{j-obj.n_int(i)})),'i_mod',j-obj.n_int(i));
                end
            end   
            name_tem=cell(obj.n_typ,1);
            for i=1:obj.n_typ
                name_tem{i}=int_info(i).TypName;
            end
            [~,id]=sort(name_tem);
            obj.mat_int=obj.mat_int(id);
            obj.n_int=obj.n_int(id);
            if ~isempty(obj.n_elm)
               obj.n_elm=obj.n_elm(id,:);
            end
%-------------------------------------------------------------------------- 
            obj.n_coord = 0;
            obj.i_coord = ones(obj.n_mod,2);
            n_tem = 0;
            for i=1:obj.n_mod
                n_tem=n_tem+obj.mod{i}.var.n_coord;
                obj.i_coord(i,2)=n_tem;
                obj.n_coord=obj.n_coord+n_tem;
                if i<obj.n_mod
                    obj.i_coord(i+1,1)=n_tem+1;
                end
            end
%--------------------------------------------------------------------------       
            for im=1:obj.n_mod
                for ip=1:numel(obj.mod{im}.prop)
                    if strcmp(obj.mod{im}.prop{ip},'needAddInit')
                        obj=obj.(['addInit_' obj.name{im}]);
                    end
                end
            end  
%--------------------------------------------------------------------------
%             obj.chemistry=chemistry(obj);
%             obj.force=force(obj);
           obj.TypForce = TypForce(obj);
           obj.TypChemistry = TypChemistry(obj);
%--------------------------------------------------------------------------          
        end
%==========================================================================        
        [f] = plot(obj,varargin);
        [obj] = addInit_ModBrownianMotor(obj,varargin);
    end
%==========================================================================
%==========================================================================  
    methods (Static)
        function [name_all,id_mod] = getModCoord(name_new,mod)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('name_new', @(x) iscell(x));
            ip.addRequired('mod', @(x) iscell(x));
            ip.parse(name_new,mod);
%--------------------------------------------------------------------------
            n_new=numel(name_new);
            n_mod=numel(mod);
            for i=1:n_new
                name_new{i}=class(name_new{i});
            end
            [name_new,~]=sort(name_new);
            id_mod=zeros(n_new,1);
            for i=1:n_new
                for j=1:n_mod
                    if strcmp(name_new{i},class(mod{j}))
                        id_mod(i)=j;
                    end
                end
            end
            name_all=name_new{1};
            if n_new>=2
            for i=2:n_new
                name_all=[name_all '_' name_new{i}];
            end
            end
%--------------------------------------------------------------------------
        end
    end
end

