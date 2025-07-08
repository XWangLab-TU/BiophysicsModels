function [f] = plot(f,m,varargin)
ip = inputParser;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addParameter('fig', [], @isobject);
ip.addParameter('FaceAlpha', 1, @isnumeric);
ip.addParameter('ScatterSize', 10, @isnumeric);
ip.addParameter('Color', [1 0 0], @isnumeric);
ip.addParameter('PaperPosition', [0 0 5 5], @isnumeric);
ip.parse(f,m,varargin{:});
%==========================================================================
% PaperPosition = ip.Results.PaperPosition;
% FaceAlpha = ip.Results.FaceAlpha;
% ScatterSize=ip.Results.ScatterSize;
% Color=ip.Results.Color;
fig = ip.Results.fig;
%==========================================================================
if isempty(fig)
    fig=figure;
else
    figure(fig); hold on;
end
% quiver3(obj.var.coord(:,1)',obj.var.coord(:,2)',obj.var.coord(:,3)',...
%         force.f_comp.ModMembrane_ModSubstrate(:,3)',...
%         force.f_comp.ModMembrane_ModSubstrate(:,4)',...
%         force.f_comp.ModMembrane_ModSubstrate(:,5)');
%f.PaperPosition = PaperPosition;

f_test=f.int_comp.ModClathrin_ModMembrane{2};
quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),f_test(:,1),f_test(:,2),f_test(:,3),'linewidth',2); 
    
    
