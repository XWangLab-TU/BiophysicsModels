function [M,breakOffInfo] = breakOff_ModClathrin(dyn,M,varargin)
%--------------------------------------------------------------------------
        % breakOff_ModMembrane performs the breakoff operation for
        % @ModClathrin during dynamics to check and add linkage if two
        % clathrin triskelia have come too close
        % input: 
        % dyn - a @dynamics object
        % M - @model object including all of the modules in dynamics
        % output:
        % breakOff - currently unset
        % optional:
        % see variable arguments
        %   See also TimeEval, breakOff_ModMembrane
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dyn', @(x) isa(x,'dynamics'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(dyn,M,varargin{:});
%--------------------------------------------------------------------------
%========================================================================================================================== 
[M.mod{M.i_mod.ModClathrin},changed] = link(M.mod{M.i_mod.ModClathrin});
breakOffInfo=struct('breakOff',changed,'id',[]);