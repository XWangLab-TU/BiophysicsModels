function [Vver] = Volume(m,varargin)
%--------------------------------------------------------------------------
        % computes the volume of vesicle given by @ModMembrane
        % input: 
        % m - a @ModMembrane object
        % output:
        % volume at each triangle: Vver
        % optional:
        % see variable arguments
        %   See also Surface
        % Author: Xinxin Wang
        % email: wangxinxin8627@gmail.com
        % date: 2025/08/05
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addParameter('idFace', [], @isnumeric);
ip.addParameter('O', [-100 -100 -100], @isnumeric);
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(m,varargin{:});
%----------------------------------------------------------------------------------------
idFace=ip.Results.idFace;
nTriAll=size(m.var.face_unq,1);
if isempty(idFace)
    nTri=nTriAll;
    idFace=1:nTri;
else
    nTri=numel(idFace);
end
O=ip.Results.O;
%----------------------------------------------------------------------------------------
iVer=m.var.face_unq(idFace,:);
A=m.var.coord(iVer(:,1),:);
B=m.var.coord(iVer(:,2),:);
C=m.var.coord(iVer(:,3),:);
Edge1=B-A; 
Edge2=B-C;
CrossRes=cross(Edge1,Edge2);
dirV=sum(CrossRes.*(A-O),2);
idtem=dirV>0;
dirV(idtem)=-1;
dirV(~idtem)=1;
Vver=zeros(nTriAll,1);
for iTri=1:nTri
    iTem=idFace(iTri);
    XYZ=m.var.coord(iVer(iTem,:),:)-O;
    Vver(iTem)=abs(-XYZ(3,1)*XYZ(2,2)*XYZ(1,3)+XYZ(2,1)*XYZ(3,2)*XYZ(1,3)+XYZ(3,1)*XYZ(1,2)*XYZ(2,3)...
                   -XYZ(1,1)*XYZ(3,2)*XYZ(2,3)-XYZ(2,1)*XYZ(1,2)*XYZ(3,3)+XYZ(1,1)*XYZ(2,2)*XYZ(3,3))/6;    
end
Vver=Vver.*dirV;
%----------------------------------------------------------------------------------------
end


