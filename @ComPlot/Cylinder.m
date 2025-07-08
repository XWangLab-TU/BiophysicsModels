function f = Cylinder(x,y,z,h,r,theta,phi, varargin)
ip = inputParser;
ip.addRequired('x', @(x) isnumeric(x));
ip.addRequired('y', @(x) isnumeric(x));
ip.addRequired('z', @(x) isnumeric(x));
ip.addRequired('h', @(x) isnumeric(x));
ip.addRequired('r', @(x) isnumeric(x));
ip.addRequired('theta', @(x) isnumeric(x));
ip.addRequired('phi', @(x) isnumeric(x));
ip.addParameter('f', [], @isobject);
ip.addParameter('nt', 20, @isnumeric);
ip.addParameter('nz', 2, @isnumeric);
ip.parse(x,y,z,h,r,theta,phi, varargin{:});
%==========================================================================
f = ip.Results.f;
if isempty(f)
    f=figure;
    figure(f); hold on;
else
    figure(f); hold on;
end
x = ip.Results.x;
y = ip.Results.y;
z = ip.Results.z;
theta = ip.Results.theta;
phi = ip.Results.phi;
n_cyl = max(size(x));
h = ip.Results.h;
r = ip.Results.r;
nz = ip.Results.nz;
nt = ip.Results.nt;
%==========================================================================
TRI = cell(n_cyl,1);
v = cell(n_cyl,1);
%==========================================================================
for i_cyl = 1:n_cyl

% Create constant vectors 
tht = linspace(0,2*pi,nt); 
l = linspace(0,h(i_cyl),nz);

% Create cylinder 
xa = repmat(r*cos(tht),nz,1); ya = repmat(r*sin(tht),nz,1); 
za = repmat(l',1,nt);

% To close the ends 
X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0]; 
Z = [za; flipud(za); za(1,:)];

% rotate
%==========================================================================
Rxt = [1 0 0; 0 cos(theta(i_cyl)) -sin(theta(i_cyl)); 0 sin(theta(i_cyl)) cos(theta(i_cyl))];
for i_col = 1:5
   xyzR = Rxt*[X(i_col,:);Y(i_col,:);Z(i_col,:)];
   X(i_col,:) = xyzR(1,:);
   Y(i_col,:) = xyzR(2,:);
   Z(i_col,:) = xyzR(3,:);
end
Rzt = [cos(phi(i_cyl)) -sin(phi(i_cyl)) 0; sin(phi(i_cyl)) cos(phi(i_cyl)) 0; 0 0 1];
for i_col = 1:5
   xyzR = Rzt*[X(i_col,:);Y(i_col,:);Z(i_col,:)];
   X(i_col,:) = xyzR(1,:);
   Y(i_col,:) = xyzR(2,:);
   Z(i_col,:) = xyzR(3,:);
end
%rotate about y axis: Ryt = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];

% translate
X = X+x(i_cyl);
Y = Y+y(i_cyl);
Z = Z+z(i_cyl);

% Draw cylinder 
[TRI{i_cyl},v{i_cyl}]= surf2patch(X,Y,Z,'triangle'); 

end

TRI_comb = [];
v_comb = [];
df = max(max(TRI{1}));
for i_cyl = 2:n_cyl
TRI{i_cyl} = TRI{i_cyl} + df*(i_cyl-1);
end

for i_cyl = 1:n_cyl
    TRI_comb = cat(1,TRI_comb,TRI{i_cyl});
    v_comb = cat(1,v_comb,v{i_cyl});
end

p = patch('Vertices',v_comb,'Faces',TRI_comb,'facecolor',[1 0.2 0.1],'facealpha',0.8); 
%view(3); grid on; axis square; title('Cylinder','FontSize',12)
p.FaceAlpha = 1;           % remove the transparency
%p.FaceColor = 'interp';    % set the face colors to be interpolated
p.LineStyle = 'none';      % remove the lines
colormap(copper)   
l = light('Position',[-0.4 0.2 0.9],'Style','infinite');
lighting gouraud
material shiny
l.Position = [-0.1 0.6 0.8];
bound = 50;
xlim([-bound bound]);
ylim([-bound bound]);
zlim([-bound bound]);
xlabel('X');
ylabel('Y');
zlabel('Z');
view(-151,30)

