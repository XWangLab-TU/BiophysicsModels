classdef ModMembrane
%==========================================================================
%--------------------------------------------------------------------------
        % ModMembrane stores the data structure of a closed membrane
        % surface and contains various functions for parameter setting,
        % variable initiation, remeshing and etc.
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
%==========================================================================    
    properties
       var %variables representing a closed membrane surface meshed into triangles
       pm %parameters
       prop %description
       failInfo %information about topological defects found
    end
%==========================================================================
%==========================================================================   
    methods
        function obj = ModMembrane(close_surf,n_ico_sphere,n_adp,varargin)
%--------------------------------------------------------------------------
        % ModMembrane is the initiation function for @ModMembrane to setup
        % parameters and variables
        % input: 
        % close_surf - currently only support 'true'
        % n_ico_sphere - index of size of initial spherical membrane, 3 is
        % used in 2022 manuscript
        % n_adp - number of adapters, currently not applied, usually = 0
        % optional:
        % see variable arguments
        %   See also SetParameter, SetVar
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('close_surf', @(x) islogical(x));
            ip.addRequired('n_ico_sphere', @(x) isnumeric(x));
            ip.addRequired('n_adp', @(x) isnumeric(x));
            ip.addParameter('unit', [], @isobject); %adjust natural unit
            ip.addParameter('update', false, @islogical);
            ip.addParameter('pm_exist', [], @isstruct);
            ip.addParameter('case_non_close', 1, @isnumeric);
            ip.parse(close_surf,n_ico_sphere,n_adp,varargin{:}); 
%--------------------------------------------------------------------------
            case_non_close=ip.Results.case_non_close;
            unit=ip.Results.unit;
            if isempty(unit)
                unit=ComUnit('erg',ComUnit.nm_to_cm(1),300,ComUnit.kBT_to_erg(1,300)); 
                warning('no unit assigned, using 1nm, 1kBT and 300K as natural units');
            end
%--------------------------------------------------------------------------
            if ip.Results.update==false
            obj.prop = {'Particle','varDt','needBreakOff','forceRelateMod'};
            obj.failInfo='';
            [obj] = SetParameter(obj,'close_surf',close_surf,'n_ico_sphere',n_ico_sphere,'unit',unit);
            if close_surf==true
                [ver_tem,face_tem] = icosphere(obj);
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
                [obj] = MeshOri(obj, obj.var.face_unq(1,:));
            else
                if case_non_close==1
                    [ver_tem,face_tem] = plane(obj);
                elseif case_non_close==2
                    [ver_tem,face_tem] = tube(obj);
                else
                    error('not the case for non-close surface');
                end
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
                [ver_tem,face_tem] = adjust_bound(obj);
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
%--------------------------------------------------------------------------
%%
                j_save=1;
                a=[obj.var.coord(obj.var.face_unq(j_save,2),1)-obj.var.coord(obj.var.face_unq(j_save,1),1),...
                    obj.var.coord(obj.var.face_unq(j_save,2),2)-obj.var.coord(obj.var.face_unq(j_save,1),2),...
                    obj.var.coord(obj.var.face_unq(j_save,2),3)-obj.var.coord(obj.var.face_unq(j_save,1),3)];
                b=[obj.var.coord(obj.var.face_unq(j_save,3),1)-obj.var.coord(obj.var.face_unq(j_save,2),1),...
                    obj.var.coord(obj.var.face_unq(j_save,3),2)-obj.var.coord(obj.var.face_unq(j_save,2),2),...
                    obj.var.coord(obj.var.face_unq(j_save,3),3)-obj.var.coord(obj.var.face_unq(j_save,2),3)];
                u_face_tem=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];
                if dot(u_face_tem,[0,0,1]) > 0
                    obj.var.face_unq(1,:)=[obj.var.face_unq(1,2),obj.var.face_unq(1,1),obj.var.face_unq(1,3)];
                end
%%                
%--------------------------------------------------------------------------
                [obj] = MeshOri(obj, obj.var.face_unq(1,:));
                obj.var.coord=obj.var.coord-mean(obj.var.coord);
            end
            else
                obj.pm=ip.Results.pm_exist;
                if isempty(ip.Results.pm_exist)
                    error('update option needs pm_exist');
                end
            if close_surf==true
                [ver_tem,face_tem] = icosphere(obj);
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
                [obj] = MeshOri(obj, obj.var.face_unq(1,:));
            else
                if case_non_close==1
                    [ver_tem,face_tem] = plane(obj);
                elseif case_non_close==2
                    [ver_tem,face_tem] = tube(obj);
                else
                    error('not the case for non-close surface');
                end
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
                [ver_tem,face_tem] = adjust_bound(obj);
                [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem,'n_adp',n_adp);
                %%
                j_save=1;
                a=[obj.var.coord(obj.var.face_unq(j_save,2),1)-obj.var.coord(obj.var.face_unq(j_save,1),1),...
                    obj.var.coord(obj.var.face_unq(j_save,2),2)-obj.var.coord(obj.var.face_unq(j_save,1),2),...
                    obj.var.coord(obj.var.face_unq(j_save,2),3)-obj.var.coord(obj.var.face_unq(j_save,1),3)];
                b=[obj.var.coord(obj.var.face_unq(j_save,3),1)-obj.var.coord(obj.var.face_unq(j_save,2),1),...
                    obj.var.coord(obj.var.face_unq(j_save,3),2)-obj.var.coord(obj.var.face_unq(j_save,2),2),...
                    obj.var.coord(obj.var.face_unq(j_save,3),3)-obj.var.coord(obj.var.face_unq(j_save,2),3)];
                u_face_tem=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];
                if dot(u_face_tem,[0,0,1]) > 0
                    obj.var.face_unq(1,:)=[obj.var.face_unq(1,2),obj.var.face_unq(1,1),obj.var.face_unq(1,3)];
                end
                [obj] = MeshOri(obj, obj.var.face_unq(1,:));
                obj.var.coord=obj.var.coord-mean(obj.var.coord);
            end
            end
        end
%==========================================================================
%==========================================================================        
%         function A = Area(obj)                     
%             pmc=zeros(10,1);
%             pmc(1) = obj.pm.nt;
%             pmc(2) = obj.pm.dt;
%             pmc(3) = obj.pm.P;
%             pmc(4) = obj.pm.k_c;
%             pmc(5) = obj.pm.k_e;
%             pmc(6) = obj.pm.dr;
%             pmc(7) = obj.pm.k_V;
%             pmc(8) = obj.pm.k_A;
%             pmc(9) = obj.pm.V0;
%             pmc(10) = obj.pm.A0;
%             pmc(11) = obj.pm.nAVmean;
%             pmc(12) = obj.pm.k_a;
%             
%             j_T = obj.var.j_T; j_T(isnan(j_T)) = 0;
%             
%             [A]=obj.AreaMex...
%                 (obj.var.coord,...
%                 pmc,...
%                 obj.var.edge_all,...
%                 obj.var.face_unq,...
%                 j_T,...
%                 obj.var.id_on_coord,...
%                 obj.var.n_node',...
%                 obj.var.T_s',...
%                 obj.var.T_e');
%         end
% %==========================================================================        
%         function V = Volume(obj)                     
%             pmc=zeros(10,1);
%             pmc(1) = obj.pm.nt;
%             pmc(2) = obj.pm.dt;
%             pmc(3) = obj.pm.P;
%             pmc(4) = obj.pm.k_c;
%             pmc(5) = obj.pm.k_e;
%             pmc(6) = obj.pm.dr;
%             pmc(7) = obj.pm.k_V;
%             pmc(8) = obj.pm.k_A;
%             pmc(9) = obj.pm.V0;
%             pmc(10) = obj.pm.A0;
%             pmc(11) = obj.pm.nAVmean;
%             pmc(12) = obj.pm.k_a;
%             
%             j_T = obj.var.j_T; j_T(isnan(j_T)) = 0;
%             
%             [V]=obj.VolumeMex...
%                 (obj.var.coord,...
%                 pmc,...
%                 obj.var.edge_all,...
%                 obj.var.face_unq,...
%                 j_T,...
%                 obj.var.id_on_face,...
%                 obj.var.n_node',...
%                 obj.var.T_s',...
%                 obj.var.T_e');
%         end        
% %==========================================================================        
%         function H = Helfrich(obj)                     
%             pmc=zeros(10,1);
%             pmc(1) = obj.pm.nt;
%             pmc(2) = obj.pm.dt;
%             pmc(3) = obj.pm.P;
%             pmc(4) = obj.pm.k_c;
%             pmc(5) = obj.pm.k_e;
%             pmc(6) = obj.pm.dr;
%             pmc(7) = obj.pm.k_V;
%             pmc(8) = obj.pm.k_A;
%             pmc(9) = obj.pm.V0;
%             pmc(10) = obj.pm.A0;
%             pmc(11) = obj.pm.nAVmean;
%             pmc(12) = obj.pm.k_a;
%             
%             j_T = obj.var.j_T; j_T(isnan(j_T)) = 0;
%             
%             [H]=obj.HelfrichMex...
%                 (obj.var.coord,...
%                 pmc,...
%                 obj.var.edge_all,...
%                 obj.var.face_unq,...
%                 j_T,...
%                 obj.var.id_on_coord,...
%                 obj.var.n_node',...
%                 obj.var.T_s',...
%                 obj.var.T_e');
%         end 
%==========================================================================        
        [f] = plot(obj,varargin);
        [ver_new,face] = plane(obj,varargin);
        [ver,face] = tube(obj,varargin);
        [ver,face] = adjust_bound(obj,varargin);
        [m] = densTran(m,k,varargin);
        [obj] = addBoundary(obj, id, varargin);
        [f] = plotFluorescence(obj,TypForce,varargin);
        [pCircle] = measureRtube(obj,n,V0,varargin);
        %[id_mesh] = getMeshID(obj,Mesh_coord, Mesh_range,Mesh_d,varargin);
        [M,remeshed] = remesh(m,M,varargin);
        [flip] = remeshFlip(m,M,idTooLong,varargin);
        [M,loc_relaxed] = remeshLocRelax(m,M,edg_add, varargin);
        [id_ring_ver1,id_ring_edg1,id_ring_ver2,id_ring_edg2] = remeshRing(m,edge_ctr, varargin);
        [split,merge] = remeshSplitMerge(m,M,id_split,id_merge,varargin);
        [m,remeshed] = remeshFlipOpt(m,j,i_edg,varargin);
        [m,remeshed,edg_add] = remeshSplitOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim,varargin);
        [m,remeshed] = remeshMergeOpt(m,j,k,id_ring_edg,edg_add_org,i_edg,rLim,varargin);
    end
%==========================================================================
%==========================================================================    
methods(Static)
    function mod_name = identify(varargin)
        mod_name = 'ModMembrane';
    end
    [mod] = docking(m1,m2,mod, varargin);
    [coord,loc_relaxed,fc,check]=locDynMex(coord,id_on_coord,A,pmc,fn,R,in,rg);
    dens=densTranMex(dens,edge_all,n_node,pm,id_on_edg);
    [var,topologicalDefect] = remeshAddVertex(pm,var,ver_new,edge_all_new,face_new, varargin);
end
end

% x=[-l*0.5,l*0.5*s3,0;...]
% l*0.5,l*0.5*s3,0;...
% -l,0,0;...
% 0,0,0;...
% l,0,0;...
% -l*0.5,-l*0.5*s3,0;...
% l*0.5,-l*0.5*s3,0;]

% f=[1 4 3;...]
% 1 2 4;...
% 2 4 5;...
% 3 4 6;...
% 6 4 7;...
% 4 5 7;]
