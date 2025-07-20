function [V] = VinMembrane2well(r1,r2,Vpm)
%--------------------------------------------------------------------------
        % VinMembrane performs the computation of the internal 
        % potential of @ModMembrane in the case of double well for flipping
        % only remeshing
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters for the internal potential
        % output:
        % V - Vin(r)
        %   See also VinMembrane2hill
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.parse(r1,r2,Vpm);
%--------------------------------------------------------------------------------------------------------
r=sqrt(sum((r1-r2).^2,2));
V0=Vpm(1);r_1=Vpm(2);r_2=Vpm(3);k_w=Vpm(4);e_b2=Vpm(5);e_w=Vpm(6);k_b21=Vpm(7);k_b22=Vpm(8);rb_21=Vpm(9);rb_22=Vpm(10);
    V = V0*(1./tanh((r-r_1)*k_w)+1./tanh((r_2-r)*k_w)...
      + e_b2/e_w./(1+exp(-k_b21*(r-rb_21)))+e_b2/e_w./(1+exp(-k_b22*(rb_22-r)))...
        );
end
