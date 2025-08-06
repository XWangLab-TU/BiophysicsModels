function [P,m] = VvolMembrane(m,Vvol,varargin)
%--------------------------------------------------------------------------
        % VvolMembrane performs the computation of the volume restriction
        % potential of @ModMembrane
        % input: 
        % m - a ModMembrane
        % Vvol.k_V - restriction strength
        % Vvol.V0 - ideal volume
        % Vedg - parameters for the internal potential
        % output:
        % P - potential energy
        % Author: Xinxin Wang
        % email: wangxinxin8627@gmail.com
        % date: 2025/08/06
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('Vvol', @(x) isstruct(x));
ip.addParameter('idFace', [], @isnumeric);
ip.parse(m,Vvol,varargin{:});
%--------------------------------------------------------------------------------------------------------
V0=Vvol.V0;
kV=Vvol.k_V;
idFace=ip.Results.idFace;
if isempty(idFace)
    [m.var.Vver] = Volume(m);
    Vver=m.var.Vver;
else
    Vver=m.var.Vver;
    [VverReplace] = Volume(m,'idFace',idFace);
    Vver(idFace)=VverReplace(idFace);
end
V=sum(Vver);
P=kV*(V-V0)^2/V0;
end