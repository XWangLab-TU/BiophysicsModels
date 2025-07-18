function [obj] = SetParameter(obj, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModFreeParticle'));
ip.addParameter('unit', [], @isobject);
ip.addParameter('n', 3, @isnumeric);
ip.addParameter('lim_xyz', [-5 5;-5 5;-5 5], @isnumeric);
ip.addParameter('D',1.0000e-06, @isnumeric); % previous (0.0001cm)^2/s diffusion coefficient
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
unit=ip.Results.unit;
%================================================================================
pm = struct(...
    'n', ip.Results.n,...
    'lim_xyz', ip.Results.lim_xyz,...
    'mu',[]...
     );
%==========================================================================
unit=convert_high_order(unit,ip.Results.D/unit.unit.kBT,2,-1); % D=mu*kBT
pm.mu=unit.unit_nat_any;
%==========================================================================            
%%
obj.pm=pm;
end