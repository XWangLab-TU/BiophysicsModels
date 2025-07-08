function [] = makeStackMov(path, varargin)
%--------------------------------------------------------------------------
        % makeStackMov performs the making of a movie in gif for a given
        % set of images
        % input: 
        % path - path to the images
        % optional:
        % see variable arguments
        %   See also savePlot, makeRotMov
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('path', @(x) ischar(x));
ip.addParameter('DelayTime', 0.5, @isnumeric);
ip.addParameter('xyzLim', [], @isnumeric);
ip.addParameter('axis_off', false, @islogical);
ip.addParameter('PaperPosition', [0.01,0.01,0.99,0.99], @isnumeric);
ip.addParameter('OuterPosition', [0,0,1,1], @isnumeric);
ip.addParameter('dt', [], @isnumeric);
ip.addParameter('view_ang',[], @isnumeric);
ip.addParameter('print_info','%.3f', @ischar);
ip.addParameter('n_file',[], @isnumeric);
ip.parse(path,varargin{:});
%--------------------------------------------------------------------------
xyzLim=ip.Results.xyzLim;
if isempty(xyzLim)
    xyzLim=[-10 10;-10 10;-10 10];
end
axis_off=ip.Results.axis_off;
PaperPosition=ip.Results.PaperPosition;
OuterPosition=ip.Results.OuterPosition;
dt=ip.Results.dt;
if isempty(dt)
    add_text=false;
else
    add_text=true;
end
view_ang=ip.Results.view_ang;
print_info=ip.Results.print_info;
n_file=ip.Results.n_file;
%--------------------------------------------------------------------------
pathFolder=[path filesep 'StackMovie'];
mkdir(pathFolder);

MyFolderInfo=dir(path);
            n_file_all=size(MyFolderInfo,1);
            if isempty(n_file)
                n_file=n_file_all;
            elseif n_file>n_file_all
                error('assigned n_file too big');              
            end
            id_fig=[];
            n_fig=0;
            for i=3:n_file
                if strcmp(MyFolderInfo(i).name(end-2:end),'fig')==true
                    n_fig=n_fig+1;
                    id_fig=[id_fig i];
                end
            end

init_gif=false;
fig_name_gif=[pathFolder filesep 'mov.gif'];
t=0;
for i=1:n_fig
    fig=openfig([path filesep MyFolderInfo(id_fig(i)).name]);
    xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
    if axis_off==true
        set(gca,'visible','off');
    end
    fig_name_png=[pathFolder filesep num2str(i) '.png'];
    if add_text
        text(xyzLim(1,2)-1,xyzLim(2,2)-1,xyzLim(3,2)-1,[num2str(t,print_info) ' s']);
        t=t+dt;
    end
    set(fig,'PaperUnits','inches','PaperPosition',PaperPosition,'OuterPosition',OuterPosition);
    if ~isempty(view_ang)
    view(view_ang);
    end
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