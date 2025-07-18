function [a]=get_a_from_ang(Phi,theta_a,phi_a)
   a=[Phi.*sin(theta_a).*cos(phi_a),Phi.*sin(theta_a).*sin(phi_a),Phi.*cos(theta_a)];
end