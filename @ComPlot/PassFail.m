function PF = PassFail(fig,path,varargin)
%--------------------------------------------------------------------------
        % PassFail determines whether a given plot matches judgement
        % input: 
        % fig - input figure object
        % path - a given path for storing result figure
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2024/01/24
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.addRequired('fig', @(x) isobject(x));
ip.addRequired('path', @(x) ischar(x));
ip.addParameter('delFigure', true, @islogical); %delete figure after judging
ip.parse(fig,path, varargin{:});
%==========================================================================
delFigure=ip.Results.delFigure;
tempFigName=['temp' num2str(int32(rand(1,1)*100000)) '.fig'];
c = uicontrol(fig,'Position',[fig.Position(1)*0.7,fig.Position(1)*0.5,fig.Position(1)*0.2,fig.Position(1)*0.05]);
c.String = 'Pass';
c.Callback = @(src,event) plotButtonPushedC(fig,path,tempFigName);

d = uicontrol(fig,'Position',[fig.Position(1)*0.7,fig.Position(1)*0.4,fig.Position(1)*0.2,fig.Position(1)*0.05]);
d.String = 'fail';
d.Callback = @(src,event) plotButtonPushedD(fig,path,tempFigName);

waitfor(fig);
fig=openfig([path, filesep tempFigName]);
Title=fig.Children(end).Title.String;
if strcmp(Title,'pass')
    PF=true;
elseif strcmp(Title,'fail')
    PF=false;
else
    error('fig title not supported');
end
close(fig);
if delFigure==true
   delete([path, filesep tempFigName]);
end
end


    function plotButtonPushedC(fig,path,tempFigName)
        title('pass');
        savefig(fig,[path, filesep tempFigName]);
    end
    function plotButtonPushedD(fig,path,tempFigName)
        title('fail');
        savefig(fig,[path, filesep tempFigName]);
    end