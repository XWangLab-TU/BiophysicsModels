function [obj] = addBoundary(obj, id, varargin)
%--------------------------------------------------------------------------
        % addBoundary performs the addition of boundary to @Modmembrane
        % input: 
        % obj - a @Modmembrane object
        % id - number of index of vertices to become the boundary
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------  
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isobject(x));
ip.addRequired('id', @(x) isnumeric(x));
ip.addParameter('update', false, @islogical);
ip.addParameter('exclusive', true, @islogical);
ip.parse(obj, id, varargin{:});
%--------------------------------------------------------------------------------------------------------
update = ip.Results.update;
exclusive=ip.Results.exclusive;
%================================================================================
%==========================================================================
Nid=numel(id);
obj.var.id_bound=id;
obj.var.n_bound=numel(obj.var.id_bound);

id_tem=zeros(size(obj.var.edge_all));
    for i=1:obj.var.n_bound
        id_tem(obj.var.edge_all==obj.var.id_bound(i))=1;
    end
    if exclusive==true
        id_tem=sum(id_tem,2)>0;
    else
        id_tem=sum(id_tem,2)==2;
    end
    obj.var.id_on_edg=(1:obj.var.n_edg)';
    obj.var.id_on_edg(id_tem)=[];
    obj.var.n_on_edg=numel(obj.var.id_on_edg);
    obj.var.id_on_coord=(1:obj.var.n_coord)';
    obj.var.id_on_coord(obj.var.id_bound)=[];
    obj.var.n_on_coord=numel(obj.var.id_on_coord);

    [ver_tem,face_tem] = adjust_bound(obj);
    [obj] = SetVar(obj,'coord',ver_tem,'face',face_tem);
    obj.var.id_bound=(1:Nid)';
    obj.var.n_bound=numel(obj.var.id_bound);

    id_tem=zeros(size(obj.var.edge_all));
    for i=1:obj.var.n_bound
        id_tem(obj.var.edge_all==obj.var.id_bound(i))=1;
    end
    if exclusive==true
        id_tem=sum(id_tem,2)>0;
    else
        id_tem=sum(id_tem,2)==2;
    end
    obj.var.id_on_edg=(1:obj.var.n_edg)';
    obj.var.id_on_edg(id_tem)=[];
    obj.var.n_on_edg=numel(obj.var.id_on_edg);
    obj.var.id_on_coord=(1:obj.var.n_coord)';
    obj.var.id_on_coord(obj.var.id_bound)=[];
    obj.var.n_on_coord=numel(obj.var.id_on_coord);
end
%==========================================================================

