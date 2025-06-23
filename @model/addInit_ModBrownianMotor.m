function [M] = addInit_ModBrownianMotor(M,varargin)
%--------------------------------------------------------------------------
        % addInit_ModBrownianMotor performs the additional initiation for
        % @ModBrownianMotor upon adding @ModBrownianMotor to @model because
        % the train and trail needs to be identified
        % input: 
        % dyn - a @dynamics object
        % M - @model object including all of the modules in dynamics
        % optional:
        % see variable arguments
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/02
%--------------------------------------------------------------------------  
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(M,varargin{:});
%--------------------------------------------------------------------------
%========================================================================================================================== 
i_mod_train=M.i_mod.(M.mod{M.i_mod.ModBrownianMotor}.var.train.modNameTrain);
i_mod_trail=M.i_mod.(M.mod{M.i_mod.ModBrownianMotor}.var.trail.modNameTrail);
idx_train=randsample(M.mod{i_mod_train}.var.n_coord,M.mod{M.i_mod.ModBrownianMotor}.var.n_coord);
idx_trail=randsample(M.mod{i_mod_trail}.var.n_coord,M.mod{M.i_mod.ModBrownianMotor}.var.n_coord);
M.mod{M.i_mod.ModBrownianMotor}=SetVar(...
M.mod{M.i_mod.ModBrownianMotor},...
M.mod{M.i_mod.ModBrownianMotor}.var.train.modNameTrain,...
M.mod{M.i_mod.ModBrownianMotor}.var.trail.modNameTrail,...
'update', true,...
'i_mod_train', i_mod_train,...
'idx_train', idx_train, ...
'i_mod_trail',i_mod_trail,...
'idx_trail', idx_trail);
M.mod{M.i_mod.ModBrownianMotor}.var.coord=M.mod{i_mod_train}.var.coord(idx_train,:);
end