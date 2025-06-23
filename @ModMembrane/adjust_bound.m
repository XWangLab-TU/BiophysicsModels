function [ver,face] = adjust_bound(obj,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isobject(x));
ip.addParameter('r', 0.15, @isnumeric);
ip.parse(obj, varargin{:});
%--------------------------------------------------------------------------
r=ip.Results.r;
%%
    max_bound=max(obj.var.id_bound);
    min_on=min(obj.var.id_on_coord);
    face=obj.var.face_unq;
    ver=obj.var.coord;
    id_bound=obj.var.id_bound;
    id_on=obj.var.id_on_coord;
    while min_on<max_bound
        face_save=face;
        ver_save=ver;
        id_bound_save=id_bound;
        id_on_save=id_on;
        face(face_save==max_bound)=min_on;
        face(face_save==min_on)=max_bound;
        ver(min_on,:)=ver_save(max_bound,:);
        ver(max_bound,:)=ver_save(min_on,:);
        id_bound(id_bound_save==max_bound)=min_on;
        id_on(id_on_save==min_on)=max_bound;
        max_bound=max(id_bound);
        min_on=min(id_on);
    end