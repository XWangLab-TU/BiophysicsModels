function [id_mesh] = getMeshID(coord,Mesh_min,Mesh_d,dim,nXYZ,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('coord', @(x) isnumeric(x));
ip.addRequired('Mesh_min', @(x) isnumeric(x));
ip.addRequired('Mesh_d', @(x) isnumeric(x));
ip.addRequired('dim', @(x) isnumeric(x));
ip.addRequired('nXYZ', @(x) isnumeric(x));
ip.parse(coord,Mesh_min,Mesh_d,dim,nXYZ,varargin{:});
%--------------------------------------------------------------------------
id_mesh=floor(coord-Mesh_min(1:dim)./Mesh_d+0.5);
ny=nXYZ(2);
nz=nXYZ(3);
id_mesh(id_mesh<0)=nan;
id_mesh=id_mesh(:,1)*ny*nz+id_mesh(:,2)*nz+id_mesh(:,3)+1;