function [dyn] = get_dt(dyn,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.parse(dyn,varargin{:});
%--------------------------------------------------------------------------
            dyn.pm.dt_var=min(dyn.pm.dt);
            dyn.pm.dt=[];