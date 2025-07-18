function [E]=getE(a,Phi)
   E=zeros(3,3,size(a,1));
   C=Phi.*sin(Phi)./(1-cos(Phi));
   factor2=(1-0.5*C)./(Phi.^2);
   E(1,1,:)=0.5*C+factor2.*a(:,1).*a(:,1);
   E(1,2,:)=0.5*a(:,3)+factor2.*a(:,1).*a(:,2);
   E(1,3,:)=-0.5*a(:,2)+factor2.*a(:,1).*a(:,3);
   E(2,1,:)=-0.5*a(:,3)+factor2.*a(:,2).*a(:,1);
   E(2,2,:)=0.5*C+factor2.*a(:,2).*a(:,2);
   E(2,3,:)=0.5*a(:,1)+factor2.*a(:,2).*a(:,3);
   E(3,1,:)=0.5*a(:,2)+factor2.*a(:,3).*a(:,1);
   E(3,2,:)=-0.5*a(:,1)+factor2.*a(:,3).*a(:,2);
   E(3,3,:)=0.5*C+factor2.*a(:,3).*a(:,3);

%    O(1,1,:)=O(1,1,:)+Phi.*cos(Phi).*factor2;
%    O(1,2,:)=O(1,2,:)-a(:,3).*sin(Phi).*factor2;
%    O(1,3,:)=O(1,3,:)+a(:,2).*sin(Phi).*factor2;
%    O(2,1,:)=O(2,1,:)+a(:,3).*sin(Phi).*factor2;
%    O(2,2,:)=O(2,2,:)+Phi.*cos(Phi).*factor2;
%    O(2,3,:)=O(2,3,:)-a(:,1).*sin(Phi).*factor2;
%    O(3,1,:)=O(3,1,:)-a(:,2).*sin(Phi).*factor2;
%    O(3,2,:)=O(3,2,:)+a(:,1).*sin(Phi).*factor2;
%    O(3,3,:)=O(3,3,:)+Phi.*cos(Phi).*factor2;
end