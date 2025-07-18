function [obj] = SetParameter(obj, varargin)
%--------------------------------------------------------------------------
        % SetParameter performs the parameter setting for @ModPolymerChain
        % allthe units are in dimensionless defined by @ComUnit
        % input: 
        % obj - a @ModPolymerChain object
        % optional:
        % see variable arguments
        %   See also SetVar
        % reference: https://doi.org/10.1016/j.procs.2017.05.152
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/03
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('obj', @(x) isa(x,'ModPolymerChain'));
ip.addParameter('unit', [], @isobject);
ip.addParameter('r', 10, @numeric); %10nm
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------------------------------------
unit=ip.Results.unit;
%================================================================================
pm=struct('r',ip.Results.r,'mu',1000,'kBT',1,'Vpm',[]);
pm.Vpm=struct('r_1',0, 'dr', 0.0001, 'r_2',2.1,... %scale
              'epsilon',1,'sigma',1.8,'LJorder',12,... %LJ potential
              'kFENE',1,'R0',2.1,...  %FENE potential
              'kappa',0,'Theta0',0.6, 'kappaDir',100,... %worm potential
              'rsq_std', 0.5, 'std_std', 0.0001 ... %stored force parameters
              );
%================================================================================
u=convert(unit,'length',ComUnit.nm_to_cm(pm.r));
pm.r=u.unit_nat.length;
%--------------------------------------------------------------------------------------------------------
%%
obj.pm=pm;
%%
end