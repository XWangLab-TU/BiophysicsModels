function [coord,dt] = translation(coord,force,mu,dt,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('coord', @(x) isnumeric(x));
            ip.addRequired('force', @(x) isnumeric(x));
            ip.addRequired('mu', @(x) isnumeric(x));
            ip.addRequired('dt', @(x) isnumeric(x));
            ip.addParameter('sMax', 0.01, @isnumeric);
            ip.parse(coord,force,mu,dt,varargin{:});
%--------------------------------------------------------------------------
sMax=ip.Results.sMax;
%--------------------------------------------------------------------------
   dt_attemp=sMax/mu;
   if dt_attemp<dt
       dt=dt_attemp;
   end
   coord=coord+force*mu*dt;                               
end