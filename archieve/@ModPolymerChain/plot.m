function [f] = plot(obj,varargin)
%--------------------------------------------------------------------------
        % plot performs the visulization of @ModPolymerChain
        % input: 
        % obj - a @ModPolymerChain object
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModPolymerChain'));
ip.addParameter('f', [], @isobject);
ip.addParameter('col', [], @isnumeric);
ip.addParameter('simple', false, @islogical);
ip.addParameter('facealpha', 1, @isnumeric);
ip.parse(obj, varargin{:}); 
%----------------------------------------------------------------------------------------
col=ip.Results.col;
if isempty(col)
    col=rand(obj.var.n_coord,3);
end
if isempty(ip.Results.f)
    f=figure;hold on;
else
    figure(ip.Results.f); hold on;
end
facealpha=ip.Results.facealpha;
%----------------------------------------------------------------------------------------
%%
if ip.Results.simple==false
[x,y,z] = sphere(10);
for i=1:obj.var.n_coord
            hold on;
            X=obj.var.coord(i,1);
            Y=obj.var.coord(i,2);
            Z=obj.var.coord(i,3);
            n = max(size(X));
            for i_plot = 1:n
                re_size=obj.pm.r;
                s1=surf(x*re_size+X(i_plot),y*re_size+Y(i_plot),z*re_size+Z(i_plot));
                %set(s1,'facecolor',col_tem,'EdgeColor','none');
                set(s1,'facecolor',col(i,:),'EdgeColor','none','facealpha',facealpha);
            end
end
else
    
end

