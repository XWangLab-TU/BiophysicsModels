function [] = makeRotMov(path, varargin)
%--------------------------------------------------------------------------
        % makeRotMov performs the making of a rotating movie in gif for a given
        % set of 3D images
        % input: 
        % path - path to the images
        % optional:
        % see variable arguments
        %   See also savePlot, makeStackMov
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('path', @(x) ischar(x));
ip.addParameter('DelayTime', 0.5, @isnumeric);
ip.parse(path,varargin{:});
%--------------------------------------------------------------------------
n_char=size(path,2);
iLastSep=0;
for i=1:n_char
    if strcmp(path(i),filesep)==1
        iLastSep=i;
    end
end
pathFolder=[path(1:iLastSep) 'RotMovie'];
mkdir(pathFolder);
n_ang=20;
d_ang=2*pi/n_ang;
view_theta=10; view_phi=(0:d_ang:2*pi-d_ang)/(2*pi)*360;

init_gif=false;
fig_name_gif=[pathFolder filesep 'mov.gif'];
for i=1:n_ang
    fig=openfig(path);
    view(view_phi(i),view_theta);
    fig_name_png=[pathFolder filesep num2str(i) '.png'];
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 10 10]);
    print(fig,fig_name_png,'-dpng','-r100');
    close(fig);
    pngImage=imread(fig_name_png);
    [gifImage,gifMap] = rgb2ind(pngImage,256);
    if init_gif == false
        imwrite(gifImage,gifMap,fig_name_gif,'gif','DelayTime',ip.Results.DelayTime, 'Loopcount',inf);
        init_gif=true;
    else
        imwrite(gifImage,gifMap,fig_name_gif,'gif','DelayTime',ip.Results.DelayTime,'WriteMode','append');
    end
end

end