function [obj] = SetParameter(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isobject(x));
ip.addParameter('unit', [], @isobject);
ip.addParameter('dt', 0.01, @isnumeric); %s^-1
ip.addParameter('konG', 32, @isnumeric); %s^-1
ip.addParameter('koff', 0.1, @isnumeric); %s^-1 
ip.addParameter('kbr', 0.2, @isnumeric); %s^-1 
ip.addParameter('ksev', 0.05, @isnumeric); %s^-1 
ip.addParameter('kcap', 0.01, @isnumeric); %s^-1 
ip.addParameter('NmaxBR', 1000, @isnumeric);
ip.addParameter('Dsu', 2.7, @isnumeric); %in nm
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
u=ip.Results.unit;
%--------------------------------------------------------------------------
%================================================================================
obj.pm = struct(...
    'dt',ip.Results.dt,...
    'konG', ip.Results.konG,...
    'koff', ip.Results.koff,...
    'kbr',ip.Results.kbr,...
    'ksev',ip.Results.ksev,...
    'kcap',ip.Results.kcap,...
    'NmaxBR', ip.Results.NmaxBR,...
    'Dsu',ip.Results.Dsu,...
    'AngBr',70/180*pi,...
    'BR_factor',[]... % su_diameter/kBT
     );
%--------------------------------------------------------------------------
u=convert(u,'length',ComUnit.nm_to_cm(obj.pm.Dsu));
obj.pm.Dsu=u.unit_nat.length;
%--------------------------------------------------------------------------
u=convert(u,'energy',ComUnit.kBT_to_erg(1,u.unit.T));
kBT=u.unit_nat.energy;
%--------------------------------------------------------------------------
obj.pm.BR_factor=obj.pm.Dsu/kBT;
end