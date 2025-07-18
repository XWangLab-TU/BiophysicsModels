function [obj] = SetVar(obj,varargin)
%--------------------------------------------------------------------------
        % SetVar performs the variable setup for @ModPolymerChain
        % input: 
        % obj - a @ModPolymerChain object
        % optional:
        % see variable arguments
        %   See also SetParameter
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/03
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModPolymerChain'));
ip.addParameter('update', false, @islogical);
ip.addParameter('nChain', [], @isnumeric);
ip.addParameter('nSubunit', [], @isnumeric);
ip.addParameter('coord1stSU', [], @isnumeric);
ip.addParameter('ang1stSU', [], @isnumeric);
ip.parse(obj,varargin{:});
%--------------------------------------------------------------------------
update = ip.Results.update;
%--------------------------------------------------------------------------
nChain=ip.Results.nChain;
nSubunit=ip.Results.nSubunit;
if isempty(nChain)
    nChain=1;
end
if isempty(nSubunit)
    nSubunit=1;
end
%--------------------------------------------------------------------------
coord1stSU=ip.Results.coord1stSU;
ang1stSU=ip.Results.ang1stSU;
if isempty(coord1stSU)
    coord1stSU=[0 0 0];
end
if isempty(ang1stSU)
    ang1stSU=[];
end
%==========================================================================
if update==false
    obj.var = struct(...
    'nChain',nChain,...
    'nSubunit',nSubunit*ones(nChain,1),...
    'coord', zeros(nChain*nSubunit(1),3),...
    'n_coord', nChain*nSubunit(1),...
    'ang', zeros(nChain*nSubunit(1),2),...
    'cosAng',nan(nChain*nSubunit(1),2),...
    'dir',nan(nChain*nSubunit(1),2),...
    'edg', []...
     );
    if strcmp(obj.kind,'dyn2')
        if nChain~=1
            error([obj.kind ' must have nChain=1']);
        end
        if (nSubunit<3) 
            error([obj.kind ' must have nSubunit>=3']);
        end
        dAng=2*pi/nSubunit(1);
        R=obj.pm.r/sin(0.5*dAng);
        for i=1:nSubunit(1)
            obj.var.ang(i,2)=0;
            obj.var.ang(i,1)=dAng*i;
            obj.var.coord(i,:)=[R*sin(obj.var.ang(i,1)),R*cos(obj.var.ang(i,1)),0];
        end
        obj.var.edg=[(1:nSubunit(1)-1)',(2:nSubunit(1))'];
    else
        error([obj.kind ' is currently not supported']);
    end
else
end
end
%==========================================================================
