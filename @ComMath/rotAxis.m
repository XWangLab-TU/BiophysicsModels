function xyzR = rotAxis(r,theta,phi, varargin)
ip = inputParser;
ip.addRequired('r', @(x) isnumeric(x));
ip.addRequired('theta', @(x) isnumeric(x));
ip.addRequired('phi', @(x) isnumeric(x));
ip.addParameter('nt', 20, @isnumeric);
ip.addParameter('nz', 2, @isnumeric);
ip.parse(r,theta,phi, varargin{:});
%==========================================================================
r = ip.Results.r;
theta = ip.Results.theta;
phi = ip.Results.phi;
%==========================================================================

% rotate theta about x axis and then phi about z
%==========================================================================
n = max(size(theta));
xyzR = r;
for i = 1:n
Rxt = [1 0 0; 0 cos(theta(i)) -sin(theta(i)); 0 sin(theta(i)) cos(theta(i))];

xyzR = Rxt*xyzR;

Rzt = [cos(phi(i)) -sin(phi(i)) 0; sin(phi(i)) cos(phi(i)) 0; 0 0 1];

xyzR = Rzt*xyzR;

end
%rotate about y axis: Ryt = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];


