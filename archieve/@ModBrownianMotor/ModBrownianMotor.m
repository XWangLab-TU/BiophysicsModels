classdef ModBrownianMotor
%==========================================================================
%--------------------------------------------------------------------------
        % ModBrownianMotor stores the data structure of Brownian Motor
        % and contains various functions for parameter setting,
        % variable initiation and etc. This object is designed to move
        % polymer chains and other molecular structures
        % currently only support flashing BM, based on:
        % The flashing Brownian ratchet and Parrondo's paradox
        % by S. N. Ethier and Jiyeon Lee
        % http://dx.doi.org/10.1098/rsos.171685
        %   See also model, ModPolymerChain
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
%==========================================================================    
    properties
       var %variables representing a Brownian Motor
       pm %parameters
       prop %description
       failInfo %information about topological defects found
       kind %give a name to differenciate motors, e.g. myosin vs dyn2
    end
%==========================================================================
%==========================================================================   
    methods
        function obj = ModBrownianMotor(typ,nMotor,modNameTrain,modNameTrail,kind,varargin)
%--------------------------------------------------------------------------
        % ModBrownianMotor is the initiation function for @ModBrownianMotor
        % to setup parameters and variables
        % input: 
        % typ - currently only support 'flash'
        % nMotor - number of motor molecule
        % modNameTrain - name of train module, e.g. 'ModPolymerChain'
        % modNameTrail - name of trail module, e.g. 'ModPolymerChain'
        % kind - name of kind of motor, e.g. 'dyn2'
        % optional:
        % see variable arguments
        %   See also SetParameter, SetVar
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/30
%--------------------------------------------------------------------------
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('typ', @(x) ischar(x));
            ip.addRequired('nMotor', @(x) isnumeric(x));
            ip.addRequired('modNameTrain', @(x) ischar(x));
            ip.addRequired('modNameTrail', @(x) ischar(x));
            ip.addRequired('kind', @(x) ischar(x));
            ip.addParameter('unit', [], @isobject); %adjust natural unit
            ip.addParameter('update', false, @islogical);
            ip.parse(typ,nMotor,modNameTrain,modNameTrail,kind,varargin{:});
%--------------------------------------------------------------------------
            unit=ip.Results.unit;   
            if isempty(unit)
                unit=ComUnit('erg',ComUnit.nm_to_cm(1),300,ComUnit.kBT_to_erg(1,300)); 
                warning('no unit assigned, using 1nm, 1kBT and 300K as natural units');
            end
%--------------------------------------------------------------------------
            obj.kind=kind;
            obj.prop = {'Particle','needAddInit'}; %'follow','needBreakOff','varDt','needAddInit'
%--------------------------------------------------------------------------
            if ip.Results.update==false
                obj.failInfo='';
                [obj] = SetParameter(obj,'unit',unit);
                [obj] = SetVar(obj,modNameTrain,modNameTrail,'n_coord',nMotor);
            end
        end     
        [obj] = SetParameter(obj,varargin);
        [obj] = SetVar(obj,varargin);
        [f] = plot(obj,varargin);
    end
%==========================================================================
%==========================================================================    
    methods(Static)
    end
end
