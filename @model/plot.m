function [f] = plot(obj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'model'));
ip.addParameter('f', [], @isobject);
ip.addParameter('col', [1 0 0], @isnumeric);
ip.addParameter('FaceAlpha', 1, @isnumeric);
ip.addParameter('LineStyle', 'none', @ischar);
ip.addParameter('LineWidth', 1, @isnumeric);
ip.addParameter('plot_mem_mesh', false, @islogical);
ip.addParameter('id_mod_mem', [], @isnumeric);
ip.parse(obj, varargin{:}); 
%----------------------------------------------------------------------------------------
id_mod_mem=obj.i_mod.ModMembrane;
%----------------------------------------------------------------------------------------
if isempty(ip.Results.f)
    f=figure;
else
    figure(ip.Results.f); hold on;
end
%----------------------------------------------------------------------------------------
if ip.Results.plot_mem_mesh
   [id_mesh] = ComMath.getMeshID(obj.mod{id_mod_mem}.var.coord(obj.mod{id_mod_mem}.var.id_on_coord,:),obj.Mesh.coord,obj.Mesh.range,obj.Mesh.d);
   scatter3(obj.Mesh.coord(id_mesh,1),obj.Mesh.coord(id_mesh,2),obj.Mesh.coord(id_mesh,3),'.');
end   
