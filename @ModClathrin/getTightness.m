function [T] = getTightness(c,M,varargin)
%--------------------------------------------------------------------------
%         getTightness performs the computation of tightness of clathrin
%         input: 
%         c - @Clathrin object
%         M - @model object
%         optional:
%         see variable arguments
%           See also getWaist, getFoot
%         Author: Xinxin Wang, Danuser Lab
%         email: wangxinxin8627@gmail.com
%         date: 2024/05/24
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(c,M,varargin{:});
%--------------------------------------------------------------------------
f=M.TypForce;
M.TypForce.pm.k_ModClathrin=[1 0];
[f,~,~] = ModClathrin(f,M);
nLeg=0;
for i=1:M.mod{M.i_mod.ModClathrin}.var.n_coord
    nLeg=nLeg+sum(sum(M.mod{M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
end
nLeg=nLeg*0.5;
x=M.mod{M.i_mod.ModClathrin}.var.coord_org;
lPL=norm(x(52,1:3)-x(8,1:3));
T=1./(sqrt(f.int_V.ModClathrin)/nLeg/lPL);