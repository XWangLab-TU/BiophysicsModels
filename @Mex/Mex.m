classdef Mex
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
 
    methods(Static)
        [f_b,f_p,u_K,kH,f_AV,V] = ModMembrane(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,dens);
        [f_c,f_m,V] = ModClathrin_ModMembrane(coord_c,coord_m,id_adp,id_rep_m_c_leg,j_T,J,pm,id_td_mem,n_node,coord_org,a);
        [f_c,V]=ModClathrin(coord,connect,pmc,coord_org,a,ctrl_pt_org);
        [f_c,f_m,V,f_adp] = ModClathrin_ModMemAdapter_ModMembrane(coord_c,coord_m,id_adp,id_rep_m_c_leg,j_T,J,pm,id_td_mem,n_node,coord_org,a);
        [f_b,f_p,u_K,kH,f_AV]=ModMemAdapter_ModMembrane(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,Cadp,Kadp);
    end
end

