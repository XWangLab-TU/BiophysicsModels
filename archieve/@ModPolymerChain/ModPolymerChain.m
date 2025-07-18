classdef ModPolymerChain
%==========================================================================
%--------------------------------------------------------------------------
        % ModPolymerChain stores the data structure of Polymer Chains
        % and contains various functions for parameter setting,
        % variable initiation and etc. This object is designed to simulate
        % polymer chains like actin, dynamin and etc.
        %   See also model,ModBrownianMotor
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
       kind %give a name to differenciate chains, e.g. actin vs dyn2
    end
%==========================================================================
%==========================================================================   
    methods
        function obj = ModPolymerChain(kind,varargin)
%--------------------------------------------------------------------------
        % ModBrownianMotor is the initiation function for @ModBrownianMotor
        % to setup parameters and variables
        % input: 
        % kind - name of kind of motor, e.g. 'actin'
        % nMotor - number of motor molecule
        % optional:
        % see variable arguments
        %   See also SetParameter, SetVar
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/30
%--------------------------------------------------------------------------
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('kind', @(x) ischar(x));
            ip.addParameter('unit', [], @isobject); %adjust natural unit
            ip.addParameter('update', false, @islogical);
            ip.addParameter('nChain', [], @isnumeric);
            ip.addParameter('nSubunit', [], @isnumeric);
            ip.parse(kind,varargin{:}); 
%--------------------------------------------------------------------------
            unit=ip.Results.unit;   
            if isempty(unit)
                unit=ComUnit('erg',ComUnit.nm_to_cm(1),300,ComUnit.kBT_to_erg(1,300)); 
                warning('no unit assigned, using 1nm, 1kBT and 300K as natural units');
            end
%-------------------------------------------------------------------------- 
            obj.prop = {'Particle'};
            obj.kind=kind;
%--------------------------------------------------------------------------
            nChain=ip.Results.nChain;
            nSubunit=ip.Results.nSubunit;
            if isempty(nChain)
                nChain=1;
            end
            if isempty(nSubunit)
                nSubunit=1;
            end            
%--------------------------------------------------------------------------
            if ip.Results.update==false
                obj.failInfo='';
                [obj] = SetParameter(obj,'unit',unit);
                [obj] = SetVar(obj,'nChain',nChain,'nSubunit',nSubunit);
            end
        end     
        [obj] = SetParameter(obj,varargin);
        [obj] = SetVar(obj,varargin);
        [f] = plot(obj,varargin);
        [obj] = getCosAng(obj,varargin);
    end
%==========================================================================
%==========================================================================    
    methods(Static)
    end
end
