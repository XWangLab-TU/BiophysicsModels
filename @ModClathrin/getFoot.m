function [foot] = getFoot(c,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addParameter('idClathrinSub', [], @isnumeric);
ip.parse(c,varargin{:});
%--------------------------------------------------------------------------
idClathrinSub=ip.Results.idClathrinSub;
if isempty(idClathrinSub)
    idClathrinSub=(1:c.var.n_coord)';
end
nSub=numel(idClathrinSub);
%--------------------------------------------------------------------------
id_tem=[49 50 51];
foot=zeros(3,3,c.var.n_coord);
for i=1:nSub
    i_c=idClathrinSub(i);
    foot(:,:,i_c)=c.get_r_from_a(c.var.coord_org(id_tem,1:3),3,c.var.O(:,:,i_c))+c.var.coord(i_c,:);
end