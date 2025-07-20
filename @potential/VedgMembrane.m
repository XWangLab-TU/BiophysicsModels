function [V] = VedgMembrane(r,Vedg)
%--------------------------------------------------------------------------
        % VinMembrane performs the computation of the internal 
        % potential of @ModMembrane
        % input: 
        % r - spatial array
        % Vedg - parameters for the internal potential
        % output:
        % V - Vin(r)
        %   See also ModMembrane8Const
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('r', @(x) isnumeric);
ip.addRequired('Vedg', @(x) isstruct);
%ip.parse(r,Vedg);
%--------------------------------------------------------------------------------------------------------
    V = Vedg.V_0*(1./tanh((r-Vedg.r_1)*Vedg.k_w)+1./tanh((Vedg.r_2-r)*Vedg.k_w)...
      - (Vedg.e_b/Vedg.e_w./(1+exp(-Vedg.k_b*(r-Vedg.rb_1)))+Vedg.e_b/Vedg.e_w./(1+exp(-Vedg.k_b*(Vedg.rb_2-r))))...
        );
end
