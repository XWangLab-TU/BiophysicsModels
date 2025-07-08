function [f,V_tot] = ModMembrane_ModSubstrate(f,mod,varargin)
%--------------------------------------------------------------------------
        % ModMembrane_ModSubstrate performs the computation of the forces
        % between @ModMembrane and @ModSubstrate
        % input: 
        % f - @TypForce
        % mod - @model object
        % output:
        % Vtot - total potential value of m
        % optional:
        % see variable arguments
        %   See also ModMembrane
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
% ip.addRequired('A', @(x) isa(x,'ModMembrane'));
% ip.addRequired('B', @(x) isobject);
% ip.addRequired('mex', @(x) isa(x,'Mex'));
%ip.addParameter('k', true, @isnumeric);
ip.parse(f,mod,varargin{:}); 
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModMembrane,mod.i_mod.ModSubstrate];  %1: membrane, 2: substrate
%----------------------------------------------------------------------------------------
int_name='ModMembrane_ModSubstrate';
%----------------------------------------------------------------------------------------
if mod.mod{i_mod(1)}.var.n_coord <= mod.mod{i_mod(2)}.var.n_coord
    i_min=1;
    n_min=mod.mod{i_mod(1)}.var.n_coord;
else
    i_min=2;
    n_min=mod.mod{i_mod(2)}.var.n_coord;
end
V_tot=0;
f.int_comp.(int_name)=cell(2,1);
f.int_comp.(int_name){1}=zeros(mod.mod{i_mod(1)}.var.n_coord,3);
f.int_comp.(int_name){2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);
switch i_min
    case 1 %fewer membrane points than substrate points
        for i_on= 1:mod.mod{i_mod(1)}.var.n_on_coord
            i=mod.mod{i_mod(1)}.var.id_on_coord(i_on);
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord-mod.mod{i_mod(1)}.var.coord(i,:));
            [~,id_tem]=min(sum(f_tem.^2,2));
            f.int_comp.(int_name){1}(i,:)=f.int_comp.(int_name){1}(i,:)+f_tem(id_tem,:);
            f.int_comp.(int_name){2}(id_tem,:)=f.int_comp.(int_name){2}(id_tem,:)-f_tem(id_tem,:);
            %--------------------------------------------------------------    
            id_tem2=mod.mod{i_mod(1)}.var.j_T(i,1:mod.mod{i_mod(1)}.var.val(i));
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord(id_tem,:)-mod.mod{i_mod(1)}.var.coord(id_tem2,:));
            f.int_comp.(int_name){1}(id_tem2,:)=f.int_comp.(int_name){1}(id_tem2,:)+f_tem;
            f.int_comp.(int_name){2}(id_tem,:)=f.int_comp.(int_name){2}(id_tem,:)+sum(f_tem);
            V_tot=V_tot+sum(0.5*f.pm.(['k_' int_name])*sqrt(sum((mod.mod{i_mod(2)}.var.coord(id_tem,:)-mod.mod{i_mod(1)}.var.coord(id_tem2,:)).^2,2)));
        end
    case 2  %fewer substrate points than membrane points
        for i= 1:mod.mod{i_mod(2)}.var.n_coord
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord(i,:)-mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.id_on_coord,:));
            [~,id_tem]=min(sum(f_tem.^2,2));
            f.int_comp.(int_name){1}(mod.mod{i_mod(1)}.var.id_on_coord(id_tem),:)=f.int_comp.(int_name){1}(mod.mod{i_mod(1)}.var.id_on_coord(id_tem),:)+f_tem(id_tem,:);
            f.int_comp.(int_name){2}(i,:)=f.int_comp.(int_name){2}(i,:)-f_tem(id_tem,:);
            %--------------------------------------------------------------
            id_tem2=mod.mod{i_mod(1)}.var.j_T(mod.mod{i_mod(1)}.var.id_on_coord(id_tem),1:mod.mod{i_mod(1)}.var.val(mod.mod{i_mod(1)}.var.id_on_coord(id_tem)));
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord(i,:)-mod.mod{i_mod(1)}.var.coord(id_tem2,:));
            f.int_comp.(int_name){1}(id_tem2,:)=f.int_comp.(int_name){1}(id_tem2,:)+f_tem;
            f.int_comp.(int_name){2}(i,:)=f.int_comp.(int_name){2}(i,:)+sum(f_tem);
            V_tot=V_tot+sum(0.5*f.pm.(['k_' int_name])*sqrt(sum((mod.mod{i_mod(2)}.var.coord(i,:)-mod.mod{i_mod(1)}.var.coord(id_tem2,:)).^2,2)));
        end
end
f.int_V.(int_name)=V_tot;
f.int_tot.(int_name)=f.int_comp.(int_name);
