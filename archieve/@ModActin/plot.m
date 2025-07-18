function [f] = plot(obj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModActin'));
ip.addParameter('f', [], @isobject);
ip.addParameter('col', [1 0 0], @isnumeric);
ip.addParameter('FaceAlpha', 1, @isnumeric);
ip.addParameter('LineStyle', 'none', @ischar);
ip.addParameter('LineWidth', 1, @isnumeric);
ip.addParameter('Simple', true, @islogical);
ip.parse(obj, varargin{:}); 
%----------------------------------------------------------------------------------------
f=ip.Results.f;
if isempty(f)
    f=figure;
else
    figure(f); hold on;
end
%----------------------------------------------------------------------------------------
if ip.Results.Simple==false
    f=ComPlot.Sphere(obj.var.coord(:,1),obj.var.coord(:,2),obj.var.coord(:,3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[1 0 0]);
else
%     iBr=1;
%     idxStart=1;
%     idxEnd=obj.var.lBr(iBr);
%     f=ComPlot.Cylinder(obj.var.coord(idxStart,1),obj.var.coord(idxStart,2),obj.var.coord(idxStart,3),...
%             (obj.var.lBr(iBr)-1)*obj.pm.Dsu,obj.pm.Dsu*0.5,obj.var.angBr(idxStart,1),obj.var.angBr(idxStart,2),'f', f);
%     for iBr=2:obj.var.Nbr
%     idxStart=idxStart+obj.var.lBr(iBr-1);
%     idxEnd=idxEnd+obj.var.lBr(iBr);
%     f=ComPlot.Cylinder(obj.var.coord(idxStart,1),obj.var.coord(idxStart,2),obj.var.coord(idxStart,3),...
%             (obj.var.lBr(iBr)-1)*obj.pm.Dsu,obj.pm.Dsu*0.5,obj.var.angBr(iBr,1),obj.var.angBr(iBr,2),'f', f); 
%     end
    idxStart=[1;cumsum(obj.var.lBr(1:end-1))+1];
    idxEnd=cumsum(obj.var.lBr);
    f=ComPlot.Sphere(obj.var.coord(idxStart,1),obj.var.coord(idxStart,2),obj.var.coord(idxStart,3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[0 1 0]);
    f=ComPlot.Sphere(obj.var.coord(idxEnd,1),obj.var.coord(idxEnd,2),obj.var.coord(idxEnd,3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[0 0 1]);
    f=ComPlot.Cylinder(obj.var.coord(idxStart,1),obj.var.coord(idxStart,2),obj.var.coord(idxStart,3),...
            (obj.var.lBr-1)*obj.pm.Dsu,obj.pm.Dsu*0.5,obj.var.angBr(:,1),obj.var.angBr(:,2),'f', f);
end


% f=ComPlot.Sphere(obj.var.coord(:,1),obj.var.coord(:,2),obj.var.coord(:,3),50,'re_size',obj.pm.su_diameter*0.5,'f', f,'col',[0 0.5 1]);

end
