function [Tau]=Tau_m(a,Phi)
   C=Phi.*sin(Phi)./(1-cos(Phi));
   Tau=(C-2)./(Phi.^2).*a;
end