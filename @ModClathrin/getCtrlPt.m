function [CtrlPt] = getCtrlPt(c,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.parse(c,varargin{:});
%--------------------------------------------------------------------------
CtrlPt=zeros(6,3,c.var.n_coord);
for i_c=1:c.var.n_coord
    CtrlPt(:,:,i_c)=c.get_r_from_a(c.var.ctrl_pt_org,6,c.var.O(:,:,i_c))+c.var.coord(i_c,:);
end