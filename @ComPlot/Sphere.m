function f = Sphere(X,Y,Z,i,varargin)
ip = inputParser;
ip.addRequired('X', @(x) isnumeric(x));
ip.addRequired('Y', @(x) isnumeric(x));
ip.addRequired('Z', @(x) isnumeric(x));
ip.addRequired('i', @(x) isnumeric(x));
ip.addParameter('f', [], @isobject);
ip.addParameter('col', [0 0 1], @isnumeric);
ip.addParameter('re_size', 1, @isnumeric);
ip.parse(X,Y,Z,i,varargin{:});
%==========================================================================
X = ip.Results.X;
Y = ip.Results.Y;
Z = ip.Results.Z;
i = ip.Results.i;
f = ip.Results.f;
if isempty(f)
    f=figure;
else
    figure(f); hold on;
end
col = ip.Results.col;
re_size = ip.Results.re_size;
n = max(size(X));
%==========================================================================
[x,y,z] = sphere(i);
figure(f);hold on;
for i_plot = 1:n
s1=surf(x*re_size+X(i_plot),y*re_size+Y(i_plot),z*re_size+Z(i_plot));
set(s1,'facecolor',col,'FaceLighting','gouraud','EdgeColor','none','facealpha',1); 
end
