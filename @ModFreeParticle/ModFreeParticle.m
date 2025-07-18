classdef ModFreeParticle
%==========================================================================
%========================================================================== 
properties
    var
    pm
    prop
    failInfo
    name
end
%==========================================================================
%==========================================================================    
methods
    function obj = ModFreeParticle(n,varargin) 
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.addRequired('n', @(x) isnumeric(x));
        ip.addParameter('unit', [], @isobject);
        ip.addParameter('lim_xyz', [-5 5;-5 5;-5 5], @isnumeric);
        ip.addParameter('coordAssigned', [], @isnumeric);
        ip.addParameter('name', [], @ischar);
        ip.parse(n,varargin{:}); 
%--------------------------------------------------------------------------
        obj.prop = {'Particle'};
        obj.failInfo='';
%--------------------------------------------------------------------------        
        obj.name=ip.Results.name;
        if isempty(obj.name)
            obj.name={'default'};
        end
%--------------------------------------------------------------------------        
        unit=ip.Results.unit;
            if isempty(unit)
                unit=ComUnit('erg',ComUnit.nm_to_cm(1),300,ComUnit.kBT_to_erg(1,300)); 
                warning('no unit assigned, using 1nm, 1kBT and 300K as natural units');
            end
%-------------------------------------------------------------------------- 
        lim_xyz=ip.Results.lim_xyz;
        coordAssigned=ip.Results.coordAssigned;
%--------------------------------------------------------------------------
        [obj] = SetParameter(obj,'n',n,'unit',unit,'lim_xyz',lim_xyz);
%--------------------------------------------------------------------------
        
        [obj] = SetVar(obj,false,'coordAssigned',coordAssigned);
    end
%==========================================================================
%==========================================================================   
        [f] = plot(obj,varargin);
end
%==========================================================================
%==========================================================================  
methods (Static)
    function mod_name = identify(varargin)
        mod_name = 'ModFreeParticle';
    end
end
%==========================================================================
%==========================================================================  
end

