function [figH] = savePlot(dir_alt, fig_name,varargin)
%--------------------------------------------------------------------------
        % savePlot performs the saving process for a given figure to a
        % given directory
        % input: 
        % dir_alt - a directory for the newly saved figure
        % fig_name - a name for the newly saved figure
        % output:
        % figH - Matlab figure handle of the newly saved figure
        % optional:
        % see variable arguments
        %   See also makeRotMov, makeStackMov
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dir_alt', @(x) ischar(x));
ip.addRequired('fig_name', @(x) ischar(x));
ip.addOptional('figH', []);
ip.addParameter('PaperPosition', [0 0 5 5], @isnumeric);
ip.addParameter('OuterPosition', [0 0 5 5], @isnumeric);

ip.parse(dir_alt, fig_name, varargin{:});


figH = ip.Results.figH;
dir_alt = ip.Results.dir_alt;
fig_name = ip.Results.fig_name;
PaperPosition = ip.Results.PaperPosition;
OuterPosition=ip.Results.OuterPosition;

if isempty(figH)
figH = findall(groot, 'Type', 'figure');
end
n_fig = max(size(figH));
%==========================================================================
for i=1:n_fig
      figure(figH(i));
      set(figH(i),'WindowStyle','normal');

      lgd=findobj(gcf,'type','legend');
      if ~isempty(lgd)
      legend('boxoff');
      legend('Location','northeast');
      end
      %figure(gcf); xlabel('');ylabel('');
      figH(i).PaperUnits = 'centimeters';
       figH(i).PaperPosition = PaperPosition;
       figH(i).OuterPosition=OuterPosition;
      %      set(gca,'fontsize',8,'Linewidth',1,'Title',[],'FontName', 'Arial');
      set(gca,'fontsize',7,'Linewidth',1,'FontName', 'Arial');
      if n_fig>1
          print(figH(i),[dir_alt filesep fig_name num2str(i)],'-dpng','-r300');
          savefig(figH(i),[dir_alt filesep fig_name num2str(i)]);
      else
          print(figH(i),[dir_alt filesep fig_name],'-dpng','-r300');
          savefig(figH(i),[dir_alt filesep fig_name]);
      end
end   

