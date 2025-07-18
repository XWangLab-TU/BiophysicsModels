function r=get_r_from_a(r,n_r,O)
   r=r';
   for i=1:n_r
       r(:,i)=O(:,:)*r(:,i);
   end
   r=r';
end