function [f] = plotFluorescence(obj,TypForce,varargin)
ip = inputParser;
ip.addRequired('obj', @(x) isobject(x));
ip.addRequired('TypForce', @(x) isobject(x));
ip.addParameter('f', [], @isobject);
% ip.addParameter('xyz', 1, @isnumeric); %1:xz projection, 2:xy
ip.addParameter('zCentr', 0, @isnumeric);
ip.addParameter('zDepth', [-1 1], @isnumeric);
ip.addParameter('xyzLim', [], @isnumeric);
ip.addParameter('dr', 0.5, @isnumeric);
ip.addParameter('kBlur', 0.1, @isnumeric);
ip.addParameter('imRange', [0 35], @isnumeric);
ip.parse(obj,TypForce, varargin{:});
%==========================================================================
vertices = obj.var.coord;
f = ip.Results.f;
if isempty(f)
    f=figure;
else
    figure(f); hold on;
end
% xyz=ip.Results.xyz;
zCentr=ip.Results.zCentr;
zDepth=ip.Results.zDepth;
kBlur=ip.Results.kBlur;
imRange=ip.Results.imRange;
%==========================================================================
xyzLim=ip.Results.xyzLim;
if isempty(xyzLim)
xyzLim=[min(vertices(:,1)), max(vertices(:,1));...
        min(vertices(:,2)), max(vertices(:,2));...
        min(vertices(:,3)), max(vertices(:,3))];
end
if xyzLim(3,1)<zCentr+zDepth(1)
    xyzLim(3,1)=zCentr+zDepth(1);
end
if xyzLim(3,2)>zCentr+zDepth(2)
    xyzLim(3,2)=zCentr+zDepth(2);
end
dr=ip.Results.dr;
nLim=[(xyzLim(1,2)-xyzLim(1,1))/dr,(xyzLim(2,2)-xyzLim(2,1))/dr,(xyzLim(3,2)-xyzLim(3,1))/dr];
nLim=floor(nLim+0.5);
[obj] = ModMembrane8Helfrich(TypForce,obj,'init',true);
%%  
F=zeros(nLim(1),nLim(2));
%     image(F);axis off;
for iz=1:nLim(3)
    zRange=[xyzLim(3,1)+(iz-1)*dr,xyzLim(3,1)+iz*dr];
    idTem=(vertices(:,3)>zRange(1)) & (vertices(:,3)<zRange(2));
    xy=floor((vertices(idTem,1:2)-xyzLim(1:2,1)')/dr+0.5);
    Axy=obj.var.f.A(idTem);
    idTem=(xy(:,1)<1) | (xy(:,1)>nLim(1)) | (xy(:,2)<1) | (xy(:,2)>nLim(2));
    xy(idTem,:)=[];
    Axy(idTem)=[];
    idx = sub2ind(size(F), xy(:,1), xy(:,2));
    F(idx)=Axy;
end
%%
nBlur=51;
nBlurHalf=floor(nBlur*0.5);
x=-nBlurHalf:nBlurHalf;
y=x';
Gblur=exp(-kBlur*(x.^2+y.^2));

Fconv = conv2(F,Gblur,'same');
% Fconv=zeros(nLim(1),nLim(2));
% for ix=1+nBlurHalf:nLim(1)-nBlurHalf
%     for iy=1+nBlurHalf:nLim(2)-nBlurHalf
%         Fconv(ix-nBlurHalf:ix+nBlurHalf,iy-nBlurHalf:iy+nBlurHalf)=Fconv(ix-nBlurHalf:ix+nBlurHalf,iy-nBlurHalf:iy+nBlurHalf)+Gblur*F(ix,iy);
%     end
% end
if max(max(Fconv)) > imRange(2)
    warning('Fconv over image range');
end
% subplot(1,2,1)
% imshow(F,[0 1.1]); 
% subplot(1,2,2)
imshow(flipud(Fconv'),imRange);
%==========================================================================
end