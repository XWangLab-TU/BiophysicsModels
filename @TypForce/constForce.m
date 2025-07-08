function [f,Vconst] = constForce(f,Mname,md,idPair,idCoord,varargin)
%--------------------------------------------------------------------------
        % constForce performs the computation of the constant force read
        % from storedForce
        % input: 
        % f - @TypForce
        % Mname - module name, i.e. 'ModMembrane'
        % md - module, i.e. @ModMembrane
        % idPair - id of which point pairs to which, i.e. edge_all in @ModMembrane
        % idCoord - id of which points to include
        % see variable arguments
        %   See also storedForce
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/17
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('Mname', @(x) ischar(x));
ip.addRequired('md', @(x) isobject(x));
ip.addRequired('idPair', @(x) isnumeric(x));
ip.addRequired('idCoord', @(x) isnumeric(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(f,Mname,md,idPair,idCoord,varargin{:});
%--------------------------------------------------------------------------------------------------------
%%
d = (md.var.coord(idPair(:,2),:) - md.var.coord(idPair(:,1),:));
r = sqrt(sum(d.^2,2));
u = d./r;
dr=f.int_stored.(Mname).rn(2)-f.int_stored.(Mname).rn(1);
i_shift=f.int_stored.(Mname).rn(1)/dr-1;
i = floor(r/dr+0.5)-i_shift;
idValid=(i<=f.int_stored.(Mname).nr);
Fconst=zeros(size(u));
Fconst(idValid,:)=f.int_stored.(Mname).fn(i(idValid)).*u(idValid,:);
f.int_const.(Mname)=zeros(md.var.n_coord,3);
Vconst=sum(f.int_stored.(Mname).Vn(i(idValid)));
for i_coord = 1:numel(idCoord)
    f.int_const.(Mname)(idCoord(i_coord),:) = f.int_const.(Mname)(idCoord(i_coord),:)-sum(Fconst(idPair(:,1)==idCoord(i_coord),:),1);
    f.int_const.(Mname)(idCoord(i_coord),:) = f.int_const.(Mname)(idCoord(i_coord),:)+sum(Fconst(idPair(:,2)==idCoord(i_coord),:),1);
end
%==============================================================================
%==============================================================================
end
%==============================================================================

