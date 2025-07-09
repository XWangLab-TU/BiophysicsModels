function [dt,da,dr] = rotationDt(obj,force,dt,kBT,varargin)
%--------------------------------------------------------------------------
        % rotationDt performs the computation of ideal time step in
        % rotation
        % input: 
        % obj - a 3D object to rotate, for example @ModClathrin
        % force - total force and torque as a ()x6 dimensional variable
        % dt - suggested time step
        % kBT - thermal temperature
        % optional:
        % see variable arguments
        %   See also TimeEval
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj', @(x) isobject(x));
            ip.addRequired('force', @(x) isnumeric(x));
            ip.addRequired('dt', @(x) isnumeric(x));
            ip.addRequired('kBT', @(x) isnumeric(x));
            ip.addParameter('idSub', [], @isnumeric);
            ip.addParameter('sMax', 0.01, @isnumeric);
            ip.parse(obj,force,dt,kBT,varargin{:});
%--------------------------------------------------------------------------
            idSub=ip.Results.idSub;
            if isempty(idSub)
                idSub=1:obj.var.n_coord;
            end
            nSub=numel(idSub);
            sMax=ip.Results.sMax;
%--------------------------------------------------------------------------
dr=zeros(nSub,6);
da=zeros(nSub,6);
         for iSub=1:nSub
             ic=idSub(iSub);
             A=obj.var.O(:,:,ic);
             E=obj.var.E(:,:,ic);
             %-------------------------------------------------------------
             mu_r=E*A*obj.pm.mu_r*A'*E';
             mu_r_s=sqrtm(mu_r);
             mu_t=E*A*obj.pm.mu_t*A'*E';
             mu_t_s=sqrtm(mu_t);
             %-------------------------------------------------------------
             Tau=force(ic,4:6);
             F=force(ic,1:3);

             Qr=(mu_r*Tau');
             Qt=(mu_t*F');

             dr(iSub,1:3)=Qt';    dr(iSub,4:6)=(mu_t_s*randn(3,1))';            
             da(iSub,1:3)=Qr';    da(iSub,4:6)=(mu_r_s*randn(3,1))';
         end
%==========================================================================        
%%
dtAttempt=[vecnorm(da(:,1:3),2,2);vecnorm(dr(:,1:3),2,2)];
dtAttempt=sMax/max(dtAttempt);
if dtAttempt<dt
    dt=dtAttempt;
end
                            
end