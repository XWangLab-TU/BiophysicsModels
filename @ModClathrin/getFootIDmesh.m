function [c] = getFootIDmesh(c,Mesh,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('Mesh', @(x) isstruct(x));
ip.parse(c,Mesh,varargin{:});
%--------------------------------------------------------------------------
coordFoot=getFoot(c);
coordFoot=coordFoot(:,:);
nTem=c.var.n_coord;
coordFoot=[reshape(coordFoot(:,1:3:end),[nTem*3,1]),reshape(coordFoot(:,2:3:end),[nTem*3,1]),reshape(coordFoot(:,3:3:end),[nTem*3,1])];

c.var.idMesh= ...
ComMath.getMeshID(coordFoot,Mesh.coord, Mesh.range,Mesh.d);