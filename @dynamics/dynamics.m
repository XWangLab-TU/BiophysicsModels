classdef dynamics
%==========================================================================
%--------------------------------------------------------------------------
        % dynamics stores the data structure and contains the functions for
        % performing the dynamics computation on a given @model object
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------  
%==========================================================================  
    properties
        pm
        i_step %current simulation step
        update
        dynTyp %0: Langevin only; 1: rotation; -1: fixed
        needFollow % a module following another module
        matchFtoMod %-1: no match found; 0: single force; >0: col number of F assigned to a mod
        updateModForce %-true: compute force also update mod
        ifVarDt
        iFon %compute only the forces needed and skip other forces
        needBreakOff %idicate whether a module needs break off function to stop within dynamics loops
        needField
        field %indicate field forces and parameters
    end
%==========================================================================
%==========================================================================    
    methods
        function obj = dynamics(n_step,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('n_step', @(x) isnumeric(x));
            ip.addParameter('cm_to_nm', 1E7, @isnumeric);
            ip.addParameter('dt_const', 0.005, @isnumeric);
            ip.parse(n_step,varargin{:});
%--------------------------------------------------------------------------
            obj.pm=struct('dt',[],'n_step',n_step,...
                          'cm_to_nm',ip.Results.cm_to_nm,...
                          'dt_var',0, ...
                          'dt_const',ip.Results.dt_const ...
                          );
            obj.i_step=1;
            obj.update=true;
            obj.field=[];
        end
%==========================================================================
        [dyn] = get_dt(dyn,varargin);
        [dyn] = Preparation(dyn,mod,Fname,varargin);
        [mod,CutOff] = TimeEval(dyn,mod,Fname,varargin);
%========================================================================== break off functions    
        [mod,breakOffInfo] = breakOff_ModMembrane(dyn,mod,varargin);
        [mod,breakOffInfo] = breakOff_ModClathrin(dyn,mod,varargin);
%==========================================================================
    end
%==========================================================================
    methods (Static)
       [coord,dt] = translation(coord,force,mu,dt,varargin);
       [obj,dt] = rotation(obj,force,dt,kBT,varargin);
       [dt,da,dr] = rotationDt(obj,force,dt,kBT,varargin);
       [M,ftot] = follow(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,varargin);
    end
end

