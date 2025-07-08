function [f] = interaction(f,mod,i_mod,int_name,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('i_mod', @(x) isnumeric(x));
ip.addRequired('int_name', @(x) ischar(x));
% ip.addRequired('A', @(x) isa(x,'ModMembrane'));
% ip.addRequired('B', @(x) isobject);
% ip.addRequired('mex', @(x) isa(x,'Mex'));
%ip.addParameter('k', true, @isnumeric);
ip.parse(f,mod,i_mod,int_name,varargin{:}); 
%----------------------------------------------------------------------------------------
if mod.mod{i_mod(1)}.var.n_coord <= mod.mod{i_mod(2)}.var.n_coord
    i_min=1;
    n_min=mod.mod{i_mod(1)}.var.n_coord;
else
    i_min=2;
    n_min=mod.mod{i_mod(2)}.var.n_coord;
end

f.int_comp.(int_name)=cell(2,1);
f.int_comp.(int_name){1}=zeros(mod.mod{i_mod(1)}.var.n_coord,3);
f.int_comp.(int_name){2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);
switch i_min
    case 1 %fewer membrane points than substrate points
        for i= 1:mod.mod{i_mod(1)}.var.n_coord
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord-mod.mod{i_mod(1)}.var.coord(i,:));
            [~,id_tem]=min(sum(f_tem.^2,2));
            f.int_comp.(int_name){1}(i,:)=f_tem(id_tem,:);
            f.int_comp.(int_name){2}(id_tem,:)=-f_tem(id_tem,:);
        end
    case 2  %fewer substrate points than membrane points
        for i= 1:mod.mod{i_mod(2)}.var.n_coord
            f_tem = f.pm.(['k_' int_name])*(mod.mod{i_mod(2)}.var.coord(i,:)-mod.mod{i_mod(1)}.var.coord);
            [~,id_tem]=min(sum(f_tem.^2,2));
            f.int_comp.(int_name){1}(id_tem,:)=f_tem(id_tem,:);
            f.int_comp.(int_name){2}(i,:)=-f_tem(id_tem,:);
        end
end