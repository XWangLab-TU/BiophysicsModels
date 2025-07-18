function [obj] = SetVar(obj,modNameTrain,modNameTrail,varargin)
%--------------------------------------------------------------------------
        % SetVar performs the variable setup for @ModBrownianMotor
        % input: 
        % obj - a @ModBrownianMotor object
        % optional:
        % see variable arguments
        %   See also SetParameter
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/30
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModBrownianMotor'));
ip.addRequired('modNameTrain', @(x) ischar(x));
ip.addRequired('modNameTrail', @(x) ischar(x));
ip.addParameter('update', false, @islogical);
ip.addParameter('n_coord', [], @isnumeric);
ip.addParameter('i_mod_train', [], @isnumeric);
ip.addParameter('idx_train', [], @isnumeric);
ip.addParameter('i_mod_trail', [], @isnumeric);
ip.addParameter('idx_trail', [], @isnumeric);
ip.parse(obj,modNameTrain,modNameTrail,varargin{:});
%--------------------------------------------------------------------------------------------------------
update = ip.Results.update;
%================================================================================
if update==false
obj.var = struct(...
    'coord', zeros(ip.Results.n_coord,3),...
    'n_coord', ip.Results.n_coord,...
    'ang', zeros(ip.Results.n_coord,2),...
    'train',[],...
    'trail',[]...
     );
obj.var.train=struct('modNameTrain',modNameTrain,'i_mod',[],'idx',[]);
obj.var.trail=struct('modNameTrail',modNameTrail,'i_mod',[],'idx',[]);
else
    obj.var.train.i_mod=ip.Results.i_mod_train;
    obj.var.train.idx=ip.Results.idx_train;
    obj.var.trail.i_mod=ip.Results.i_mod_trail;
    obj.var.trail.idx=ip.Results.idx_trail;
end
end
%==========================================================================
