function [f] = plot(obj,varargin)
%--------------------------------------------------------------------------
        % plot performs the visulization of @ModFreeParticle
        % input: 
        % obj - a @ModFreeParticle object
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.addRequired('obj', @(x) isa(x,'ModFreeParticle'));
ip.addParameter('f', [], @isobject);
ip.addParameter('FaceAlpha', 1, @isnumeric);
ip.addParameter('ScatterSize', 20, @isnumeric);
ip.addParameter('Color', [1 0 0], @isnumeric);
ip.addParameter('PaperPosition', [0 0 5 5], @isnumeric);
ip.parse(obj, varargin{:});
%==========================================================================
PaperPosition = ip.Results.PaperPosition;
FaceAlpha = ip.Results.FaceAlpha;
ScatterSize=ip.Results.ScatterSize;
Color=ip.Results.Color;
f = ip.Results.f;
%==========================================================================
if isempty(f)
    f=figure;
else
    figure(f); hold on;
end
scatter3(obj.var.coord(:,1),obj.var.coord(:,2),obj.var.coord(:,3),ScatterSize,'filled','MarkerFaceColor',Color,'MarkerFaceAlpha',FaceAlpha);
f.PaperPosition = PaperPosition;