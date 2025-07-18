function [] = TestCMEforce(mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('nS', 200, @isnumeric);
ip.parse(mod,varargin{:});
%--------------------------------------------------------------------------
nS=ip.Results.nS;
dyn=dynamics(nS);
Fname={'ModClathrin';'ModMembrane';'ModMembrane_ModSubstrate';'ModClathrin_ModFreeParticle';'ModFreeParticle_ModMembrane'};
[dyn] = Preparation(dyn,mod,Fname);
%--------------------------------------------------------------------------
mod.TypForce.pm.k_ModClathrin(1)=10.; %20
mod.TypForce.pm.k_ModClathrin(2)=1000; %200
mod.TypForce.pm.k_ModClathrin_ModFreeParticle=0;
mod.TypForce.pm.k_ModFreeParticle_ModMembrane=0;
for im=1:mod.n_mod
    mod.mod{mod.i_mod.(mod.name{im})}.pm.mu=1e-30;
end
mod.mod{mod.i_mod.ModClathrin}.pm.mu=1000;
%--------------------------------------------------------------------------
[~] = TimeEval(dyn,mod,Fname,'plot_or_not', true);
%--------------------------------------------------------------------------
mod.TypForce.pm.k_ModClathrin(1)=0.; %20
mod.TypForce.pm.k_ModClathrin(2)=0; %200
mod.TypForce.pm.k_ModClathrin_ModFreeParticle=100;
mod.TypForce.pm.k_ModFreeParticle_ModMembrane=0;
for im=1:mod.n_mod
    mod.mod{mod.i_mod.(mod.name{im})}.pm.mu=1e-30;
end
mod.mod{mod.i_mod.ModFreeParticle}.pm.mu=1000000;
%--------------------------------------------------------------------------
[~] = TimeEval(dyn,mod,Fname,'plot_or_not', true,'sMax', 0.01);
%--------------------------------------------------------------------------
mod.TypForce.pm.k_ModClathrin(1)=0.; %20
mod.TypForce.pm.k_ModClathrin(2)=0; %200
mod.TypForce.pm.k_ModClathrin_ModFreeParticle=0;
mod.TypForce.pm.k_ModFreeParticle_ModMembrane=1000;
for im=1:mod.n_mod
    mod.mod{mod.i_mod.(mod.name{im})}.pm.mu=1e-30;
end
mod.mod{mod.i_mod.ModFreeParticle}.pm.mu=1000000;
%--------------------------------------------------------------------------
[~] = TimeEval(dyn,mod,Fname,'plot_or_not', true,'sMax', 0.01);
%--------------------------------------------------------------------------
mod.TypForce.pm.k_ModClathrin(1)=10.; %20
mod.TypForce.pm.k_ModClathrin(2)=1000; %200
mod.TypForce.pm.k_ModClathrin_ModFreeParticle=1000;
mod.TypForce.pm.k_ModFreeParticle_ModMembrane=1000;
% for im=1:mod.n_mod
%     mod.mod{mod.i_mod.(mod.name{im})}.pm.mu=1e-30;
% end
mod.mod{mod.i_mod.ModClathrin}.pm.mu=1000;
mod.mod{mod.i_mod.ModFreeParticle}.pm.mu=1000000;
mod.mod{mod.i_mod.ModMembrane}.pm.mu=10000;
%--------------------------------------------------------------------------
[~] = TimeEval(dyn,mod,Fname,'plot_or_not', true,'sMax', 0.01);
%--------------------------------------------------------------------------