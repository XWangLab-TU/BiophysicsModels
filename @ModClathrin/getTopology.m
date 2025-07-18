function [topology,n_ring] = getTopology(c,junction,varargin)
%--------------------------------------------------------------------------
        % init performs the topological evaluation of @ModClathrin, getting
        % whether a clathrin could close a ring 
        % input: 
        % c - @ModClathrin object
        % junction - clathrin and leg number of a pair of clathrin
        % output:
        % topology - close ring or not
        % n_ring - number of clathrin in 2 closed rings, left or right hand
        % optional:
        % see variable arguments
        %   See also add, init
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('junction', @(x) isnumeric(x));
ip.parse(c,junction,varargin{:});
%--------------------------------------------------------------------------
%%
i_c=junction(1,1);
    i_leg=junction(1,2);   
%%
% i_c_try=i_c_check_left; i_leg_try=i_leg_check_left;
% i_c_try=i_c_check_right; i_leg_try=i_leg_check_right;
% i_c_try=i_c_check_right_alt; i_leg_try=i_leg_check_right_alt;
% i_c_try=i_c_check_new; i_leg_try=i_leg_check_new;
% i_c_try=junction(2,1); i_leg_try=junction(2,2);
% fig=figure('units','normalized','outerposition',[0 0 1 1]); 
% plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',true,...
%      'simple',true,'iC_iLeg',[i_c i_leg;i_c_try i_leg_try]);%'iC_iLeg',[6 1;4 2]
%%
         topology_found=false;
         topology=false;
         n_ring=0;
         if size(junction,1)==2
             %%
             i_c_check_left=c.var.connect(i_leg,1,i_c);
             i_leg_check_left=c.var.connect(i_leg,2,i_c);
             i_leg_check_left=left_hand(i_leg_check_left);
             n_ring_left=i_c_check_left;
             
             i_c_check_left_alt=i_c_check_left;
             i_leg_check_left_alt=i_leg_check_left;
             
             i_c_check_right=c.var.connect(i_leg,1,i_c);
             i_leg_check_right=c.var.connect(i_leg,2,i_c);
             i_leg_check_right=right_hand(i_leg_check_right);
             
             i_c_check_right_alt=i_c_check_right;
             i_leg_check_right_alt=i_leg_check_right;
             n_ring_right=i_c_check_right;
             %%
             while (topology_found==false) 
              %%
                if i_c_check_left>0
                i_c_check_new=c.var.connect(i_leg_check_left,1,i_c_check_left);
                i_leg_check_new=c.var.connect(i_leg_check_left,2,i_c_check_left);
                i_leg_check_new=left_hand(i_leg_check_new);
                
                i_c_check_left=i_c_check_new;
                i_leg_check_left=i_leg_check_new;
                n_ring_left=[n_ring_left;i_c_check_left];
                end
                %%
                if ((i_c_check_left==junction(2,1)) && (i_leg_check_left==junction(2,2))) || ((i_c_check_left==junction(1,1)) && (i_leg_check_left==junction(1,2)))
                    topology_found=true;
                    topology=true;
                    n_ring=numel(unique(n_ring_left));
                end
                if i_c_check_left_alt>0
                i_c_check_new=c.var.connect(i_leg_check_left_alt,1,i_c_check_left_alt);
                i_leg_check_new=c.var.connect(i_leg_check_left_alt,2,i_c_check_left_alt);
                i_leg_check_new=right_hand(i_leg_check_new);
                
                i_c_check_left_alt=i_c_check_new;
                i_leg_check_left_alt=i_leg_check_new;
                end
                if ((i_c_check_left_alt==junction(2,1)) && (i_leg_check_left_alt==junction(2,2))) || ((i_c_check_left_alt==junction(1,1)) && (i_leg_check_left_alt==junction(1,2)))
                    topology_found=true;
                    topology=false;
                end
                %%
                if i_c_check_right>0
                i_c_check_new=c.var.connect(i_leg_check_right,1,i_c_check_right);
                i_leg_check_new=c.var.connect(i_leg_check_right,2,i_c_check_right);
                i_leg_check_new=right_hand(i_leg_check_new);
                
                i_c_check_right=i_c_check_new;
                i_leg_check_right=i_leg_check_new;
                n_ring_right=[n_ring_right;i_c_check_right];
                end
                %%
                if ((i_c_check_right==junction(2,1)) && (i_leg_check_right==junction(2,2))) || ((i_c_check_right==junction(1,1)) && (i_leg_check_right==junction(1,2)))
                    topology_found=true;
                    topology=true;
                    n_ring=numel(unique(n_ring_right));
                end
                if i_c_check_right_alt>0
                i_c_check_new=c.var.connect(i_leg_check_right_alt,1,i_c_check_right_alt);
                i_leg_check_new=c.var.connect(i_leg_check_right_alt,2,i_c_check_right_alt);
                i_leg_check_new=left_hand(i_leg_check_new);
                
                i_c_check_right_alt=i_c_check_new;
                i_leg_check_right_alt=i_leg_check_new;
                end
                if ((i_c_check_right_alt==junction(2,1)) && (i_leg_check_right_alt==junction(2,2))) || ((i_c_check_right_alt==junction(1,1)) && (i_leg_check_right_alt==junction(1,2)))
                    topology_found=true;
                    topology=false;
                end
            end
         elseif size(junction,1)==1
             %%
             topology=false;
             n_ring=[1;1];
%-------------------------------------------------------------------------- left             
             check_terminate=false;
             i_leg_check_new=left_hand(i_leg);
             i_c_check_left=c.var.connect(i_leg_check_new,1,i_c);
             i_leg_check_left=c.var.connect(i_leg_check_new,2,i_c);
             if c.var.connect(i_leg_check_new,1,i_c)>0
                 n_ring(1)=n_ring(1)+1;
             else
                 check_terminate=true;
             end
%%
             while (check_terminate==false) 
                 %%
               i_leg_check_new=left_hand(i_leg_check_left);
               i_c_check_new=c.var.connect(i_leg_check_new,1,i_c_check_left);
               i_leg_check_new=c.var.connect(i_leg_check_new,2,i_c_check_left);
               
               i_c_check_left=i_c_check_new;
               i_leg_check_left=i_leg_check_new;
               
               i_leg_switch=left_hand(i_leg_check_left);

               if i_c_check_left>0
                   n_ring(1)=n_ring(1)+1;
                   if c.var.connect(i_leg_switch,1,i_c_check_left)==0
                      check_terminate=true;
                   end
               else
                   check_terminate=true;
               end
               %%
             end
%-------------------------------------------------------------------------- right
%%
             check_terminate=false;
             i_leg_check_new=right_hand(i_leg);
             i_c_check_right=c.var.connect(i_leg_check_new,1,i_c);
             i_leg_check_right=c.var.connect(i_leg_check_new,2,i_c);
             if c.var.connect(i_leg_check_new,1,i_c)>0
                 n_ring(2)=n_ring(2)+1;
             else
                 check_terminate=true;
             end
%%
             while (check_terminate==false) 
                 %%
               i_leg_check_new=right_hand(i_leg_check_right);
               i_c_check_new=c.var.connect(i_leg_check_new,1,i_c_check_right);
               i_leg_check_new=c.var.connect(i_leg_check_new,2,i_c_check_right);
               
               i_c_check_right=i_c_check_new;
               i_leg_check_right=i_leg_check_new;
               
               i_leg_switch=right_hand(i_leg_check_right);
               if i_c_check_right>0
                   n_ring(2)=n_ring(2)+1;
                   if c.var.connect(i_leg_switch,1,i_c_check_right)==0
                      check_terminate=true;
                   end
               else
                   check_terminate=true;
               end
             end
         else
             error('junction is wrong');
         end
            
%%
%fig=figure;
% % plot(c,'f',fig,'proxy_only',false,'col',col,'simple',true,'iC_iLeg',[junction(1,1) junction(1,2)]);%'iC_iLeg',[6 1;4 2]
% subplot(1,2,1)
% plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',true,'col',col,'simple',true,'iC_iLeg',[5 1;i_c_check_left i_leg_switch]);
% subplot(1,2,2)
% plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',true,'col',col,'simple',true,'iC_iLeg',[5 1;i_c_check_left i_leg_check_left]);
% subplot(1,2,1)
% plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',true,'col',col,'simple',true,'iC_iLeg',[5 1;i_c_check_right i_leg_switch]);
% subplot(1,2,2)
% plot(mod.mod{mod.i_mod.ModClathrin},'f',fig,'proxy_only',true,'col',col,'simple',true,'iC_iLeg',[5 1;i_c_check_right i_leg_check_right]);
end

function [j]=right_hand(i)
   j=i+1;
   if j>3
       j=1;
   end
end
function [j]=left_hand(i)
   j=i-1;
   if j<1
       j=3;
   end
end