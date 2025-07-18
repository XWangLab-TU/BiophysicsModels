function [obj] = SetVar(obj,update,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModFreeParticle'));
ip.addRequired('update', @(x) islogical(x));
ip.addParameter('coordAssigned', [], @isnumeric);
ip.addParameter('coord_org', [], @isnumeric);
ip.parse(obj,update,varargin{:});
%----------------------------------------------------------------------------------------
obj.var=struct('coord',zeros(obj.pm.n,3),'n_coord',obj.pm.n,'idMesh',zeros(obj.pm.n,1),'follow',[]);
%----------------------------------------------------------------------------------------    
coordAssigned=ip.Results.coordAssigned;
if isempty(coordAssigned)
    obj.var.coord=[obj.pm.lim_xyz(1,1)+(rand(obj.pm.n,1))*(obj.pm.lim_xyz(1,2)-obj.pm.lim_xyz(1,1)),...
                   obj.pm.lim_xyz(2,1)+(rand(obj.pm.n,1))*(obj.pm.lim_xyz(2,2)-obj.pm.lim_xyz(2,1)),...
                   obj.pm.lim_xyz(3,1)+(rand(obj.pm.n,1))*(obj.pm.lim_xyz(3,2)-obj.pm.lim_xyz(3,1))];
else
    obj.var.coord=coordAssigned;
    obj.var.n_coord=size(obj.var.coord,1);
end
%----------------------------------------------------------------------------------------
obj.var.follow=struct('nameV','ModClathrin_ModFreeParticle','nameIdNeighbor','j_T','idModFollower',2,'idModFollowee',3,...
                      'idCoordFollowee',[]);
%----------------------------------------------------------------------------------------
end

