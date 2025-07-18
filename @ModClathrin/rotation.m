function [coord_final]=rotation(coord,theta,phi)
            coord = coord';
            n_tem = size(coord,2);
            n = max(size(theta));
            coord_final = zeros(3,n_tem);
            for i_tem = 1:n_tem
            xyzR = coord(:,i_tem);          
            for i = 1:n
                Rxt = [1 0 0; 0 cos(theta(i)) -sin(theta(i)); 0 sin(theta(i)) cos(theta(i))];
                
                xyzR = Rxt*xyzR;
                
                Rzt = [cos(phi(i)) -sin(phi(i)) 0; sin(phi(i)) cos(phi(i)) 0; 0 0 1];
                
                xyzR = Rzt*xyzR;
                
            end
            coord_final(:,i_tem) = xyzR;
            end
            coord_final=coord_final';
        end