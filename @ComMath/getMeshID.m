function [id_mesh] = getMeshID(coord,Mesh_coord, Mesh_range,Mesh_d,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('coord', @(x) isnumeric(x));
ip.addRequired('Mesh_coord', @(x) isnumeric(x));
ip.addRequired('Mesh_range', @(x) iscell(x));
ip.addRequired('Mesh_d', @(x) isnumeric(x));
ip.parse(coord,Mesh_coord, Mesh_range,Mesh_d,varargin{:});
%--------------------------------------------------------------------------
id_mesh=floor((coord-[Mesh_range{1}(1),Mesh_range{2}(1),Mesh_range{3}(1)])./Mesh_d+0.5);
ny=numel(Mesh_range{2});
nz=numel(Mesh_range{3});
id_mesh(id_mesh<0)=nan;
id_mesh=id_mesh(:,1)*ny*nz+id_mesh(:,2)*nz+id_mesh(:,3)+1;