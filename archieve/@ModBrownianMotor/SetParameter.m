function [obj] = SetParameter(obj, varargin)
%--------------------------------------------------------------------------
        % SetParameter performs the parameter setting for @ModBrownianMotor
        % allthe units are in dimensionless defined by @ComUnit
        % input: 
        % obj - a @ModBrownianMotor object
        % optional:
        % see variable arguments
        %   See also SetVar
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/30
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('obj', @(x) isa(x,'ModBrownianMotor'));
ip.addParameter('unit', [], @isobject);
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
unit=ip.Results.unit;
%================================================================================
pm=struct('A',0.1,'Adrag',1,'mu',1000);
%--------------------------------------------------------------------------------------------------------
%%
obj.pm=pm;
%%
end