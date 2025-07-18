function [O]=Omega(a,Phi)
   O=zeros(3,3,size(a,1));
   factor1=(1-cos(Phi))./Phi.^2;
   factor2=1./Phi;
   O(1,1,:)=a(:,1).*a(:,1).*factor1+Phi.*cos(Phi).*factor2;
   O(1,2,:)=a(:,1).*a(:,2).*factor1-a(:,3).*sin(Phi).*factor2;
   O(1,3,:)=a(:,1).*a(:,3).*factor1+a(:,2).*sin(Phi).*factor2;
   O(2,1,:)=a(:,2).*a(:,1).*factor1+a(:,3).*sin(Phi).*factor2;
   O(2,2,:)=a(:,2).*a(:,2).*factor1+Phi.*cos(Phi).*factor2;
   O(2,3,:)=a(:,2).*a(:,3).*factor1-a(:,1).*sin(Phi).*factor2;
   O(3,1,:)=a(:,3).*a(:,1).*factor1-a(:,2).*sin(Phi).*factor2;
   O(3,2,:)=a(:,3).*a(:,2).*factor1+a(:,1).*sin(Phi).*factor2;
   O(3,3,:)=a(:,3).*a(:,3).*factor1+Phi.*cos(Phi).*factor2;

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