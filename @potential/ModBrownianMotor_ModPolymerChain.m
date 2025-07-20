function [V] = ModBrownianMotor_ModPolymerChain(bm,pc,varargin)
%--------------------------------------------------------------------------
        % Vmotor performs the computation of motor potential between a
        % motor molecule (positioned at the motor's train) and a polymer chain 
        % (positioned at the motor's trail)
        % input: 
        % bm - @ModBrownianMotor
        % pc - @ModPolymerChain 
        % output:
        % V - potential array
        % optional:
        % see variable arguments
        %   See also Vin
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/03
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('bm', @(x) isa(x,'ModBrownianMotor'));
ip.addRequired('pc', @(x) isa(x,'ModPolymerChain'));
ip.addParameter('plot_or_not', false, @islogical); %whether to plot result
ip.parse(bm,pc,varargin{:});
%----------------------------------------------------------------------------------------
plot_or_not=ip.Results.plot_or_not;
%----------------------------------------------------------------------------------------
%% computation

%--------------------------------------------------------------------------    
end
