function [V] = VinMembrane2hill(r1,r2,Vpm)
%--------------------------------------------------------------------------
        % VinMembrane2hill performs the computation of the internal 
        % potential of @ModMembrane in the case of double hill for 
        % splitting and merging based remeshing
        % input: 
        % r1 - position of 1st point
        % r2 - position of 2nd point
        % Vpm - parameters for the internal potential
        % output:
        % V - Vin(r)
        %   See also VinMembrane2well
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
V0=Vpm(1);
r_1=Vpm(2);r_2=Vpm(3);
k_w=Vpm(4);
e_b1=Vpm(5);e_b2=Vpm(6);
e_w=Vpm(7);
k_b11=Vpm(8);k_b12=Vpm(9);k_b21=Vpm(10);k_b22=Vpm(11);
rb_11=Vpm(12);rb_12=Vpm(13);rb_21=Vpm(14);rb_22=Vpm(15);
    V = V0*(1./tanh((r-r_1)*k_w)+1./tanh((r_2-r)*k_w)...
      + e_b1/e_w./(1+exp(-k_b11*(r-rb_11)))+e_b1/e_w./(1+exp(-k_b12*(rb_12-r)))...
      + e_b2/e_w./(1+exp(-k_b21*(r-rb_21)))+e_b2/e_w./(1+exp(-k_b22*(rb_22-r)))...
        );
   V=V-min(V); 
end
