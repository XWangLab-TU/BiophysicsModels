function [var,topologicalDefect] = remeshAddVertex(pm,var,ver_new,edge_all_new,face_new, varargin)
%--------------------------------------------------------------------------
        % remesh performs the vertex adjustment after e.g. remeshing
        % input: 
        % pm - parameter from @ModMembrane
        % var - variable from @ModMembrane
        % ver_new - updated variable
        % edge_all_new - updated edge_all variable from @ModMembrane
        % face_new - updated face_unq variable from @ModMembrane
        % output:
        % topologicalDefect - true/false whether topologicalDefect detected
        % optional:
        % see variable arguments
        %   See also remesh
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('pm', @(x) isstruct(x));
ip.addRequired('var', @(x) isstruct(x));
ip.addRequired('ver_new', @(x) isnumeric(x));
ip.addRequired('edge_all_new', @(x) isnumeric(x));
ip.addRequired('face_new', @(x) isnumeric(x));
ip.addParameter('var_old', [], @isstruct);
ip.addParameter('dens', [], @isnumeric);
ip.parse(pm, var,ver_new, edge_all_new,face_new,varargin{:});
%--------------------------------------------------------------------------------------------------------
%================================================================================
    %%
     n_old=size(var.coord,1);
     var.coord = ver_new; 
     n_new=size(var.coord,1);
     dens=ip.Results.dens;
     if isempty(dens)
         dens=ones(n_new,1);
     end
     var.dens=dens;
     var.n_coord = n_new;
     var.face_unq = face_new;
     var.face = nan(pm.n_val_max,3,var.n_coord);
     var.val = zeros(var.n_coord,1);
     for i = 1:var.n_coord
         id_tem = var.face_unq == i;
         id_tem = sum(id_tem,2);
         var.val(i) = sum(id_tem);
         var.face(1:var.val(i),:,i) = face_new(logical(id_tem),:);
     end
     
     var.edge_all=edge_all_new;
     var.n_edg = size(var.edge_all,1);
     
%--------------------------------------------------------------------------
if pm.close_surf==true
var.id_on_edg=(1:var.n_edg)';
var.id_on_coord=(1:var.n_coord)';
var.n_on_coord=numel(var.id_on_coord);
else
    id_tem=zeros(size(var.edge_all));
    for i=1:var.n_bound
        id_tem(var.edge_all==var.id_bound(i))=1;
    end
    id_tem=sum(id_tem,2)==2;
    var.id_on_edg=(1:var.n_edg)';
    var.id_on_edg(id_tem)=[];
    var.n_on_edg=numel(var.id_on_edg);
    var.id_on_coord=(1:var.n_coord)';
    var.id_on_coord(var.id_bound)=[];
    var.n_on_coord=numel(var.id_on_coord);
end
% patch('Vertices',var.coord,'Faces',var.face_unq,'FaceVertexCData',rand(var.n_coord,3),'FaceColor','interp','facealpha',1);
%--------------------------------------------------------------------------
     %%
     [var.j_T,topologicalDefect] = comp_j_T(var,pm);
     if topologicalDefect==false
     %%
     var.j_idx = cell(var.n_coord,1);
     var.faces_T = cell(var.n_coord,1);
     for i = 1:var.n_coord
         id_temp = (var.face_unq(:,:) == i) ;
         id_temp = sum(id_temp,2) > 0;
         var.faces_T{i} = var.face_unq(id_temp,:);
         var.j_idx{i} = 1:var.val(i)+1;
         var.j_idx{i}(end) = 1;
     end
     var.J = [];
     var.T_s = size(var.n_coord,1);
     var.T_e = size(var.n_coord,1);
     for i = 1:var.n_coord
         id_temp = ~isnan(var.j_T(i,:));
         var.T_s(i) = size(var.J,2)+1;
         var.J = [var.J,var.j_T(i,id_temp)];
         var.T_e(i) = size(var.J,2);
     end
     var.J = var.J';
     var.n_node = var.T_e-var.T_s+1;
     var.u1 = cell(var.n_coord,1);var.u2 = cell(var.n_coord,1);
     for i = 1:var.n_coord
         var.u1{i} = eye(var.n_node(i));
         var.u2{i} = circshift(var.u1{i},var.n_node(i)-1);
     end
     else
     end
%==========================================================================
end
%==========================================================================

function [j_T,topologicalDefect] = comp_j_T(var,pm)
j_T = nan(var.n_coord,pm.n_val_max);
topologicalDefect=false;
for i = 1:var.n_coord
    x_tem=var.face(:,:,i);
    x_tem=x_tem(x_tem~=0);
    ic = mode(x_tem,'all');
    face_tem = var.face(1,:,i);
    face_tem = face_tem(face_tem ~= ic);
    j_T(i,1) = face_tem(1);
    j_T(i,var.val(i)) = face_tem(2);
    face_tem = var.face(:,:,i);
    for j = 2:var.val(i)-1
        j_pre = j-2; if j_pre < 1; j_pre = j_pre+var.val(i); end
        id_tem = sum((face_tem == j_T(i,j-1)),2)+sum((face_tem ~= j_T(i,j_pre)),2);
        j_T_tem = face_tem(id_tem==4,:);
        if ~isempty(j_T_tem)
            %check for topological defect
            j_T_tem=j_T_tem((j_T_tem~=ic)&(j_T_tem~=j_T(i,j-1)));
            if numel(j_T_tem)==1
                j_T(i,j) = j_T_tem;
            else
                topologicalDefect=true;
                warning('topological defect detected');
            end
        end
        if topologicalDefect==true
            break;
        end
    end
    if topologicalDefect==true
       break;
    end
end
if topologicalDefect==true
   j_T = nan(var.n_coord,pm.n_val_max);
end
end






