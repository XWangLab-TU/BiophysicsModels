function [obj] = SetParameter(obj, varargin)
%--------------------------------------------------------------------------
        % SetParameter sets parameters for @TypForce
        % input: 
        % obj - @TypForce object
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'TypForce'));
ip.addParameter('unit', [], @isobject);
ip.addParameter('k_ModMembrane_ModSubstrate', 1, @isnumeric);
ip.addParameter('k_ModClathrin_ModMembrane', 10, @isnumeric);
ip.addParameter('k_ModClathrin', [200 2000], @isnumeric);
ip.addParameter('k_ModClathrin_ModMemAdapter_ModMembrane', 10, @isnumeric);
ip.addParameter('k_ModClathrin_ModFreeParticle', 10, @isnumeric);
ip.addParameter('k_ModFreeParticle_ModMembrane', [10 0.5], @isnumeric);
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
unit=ip.Results.unit;
%================================================================================
unit=convert_high_order(unit,unit.unit.kBT/unit.unit.length^2,-2,1);
            k=unit.unit_nat_any;

pm = struct(...
    'k_ModMembrane_ModSubstrate', ip.Results.k_ModMembrane_ModSubstrate*k,...
    'k_ModClathrin_ModMembrane', ip.Results.k_ModClathrin_ModMembrane*k,...
    'k_ModClathrin', ip.Results.k_ModClathrin*k,...
    'k_ModClathrin_ModMemAdapter_ModMembrane',ip.Results.k_ModClathrin_ModMemAdapter_ModMembrane*k,...
    'k_ModClathrin_ModFreeParticle', ip.Results.k_ModClathrin_ModFreeParticle*k,...
    'k_ModFreeParticle_ModMembrane', ip.Results.k_ModFreeParticle_ModMembrane*k...
     );
obj.pm=pm;