function [waist] = getWaist(c,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.parse(c,varargin{:});
%--------------------------------------------------------------------------
id_tem=[8 16 24];
waist=zeros(3,3,c.var.n_coord);
for i_c=1:c.var.n_coord
    waist(:,:,i_c)=c.get_r_from_a(c.var.coord_org(id_tem,1:3),3,c.var.O(:,:,i_c))+c.var.coord(i_c,:);
end