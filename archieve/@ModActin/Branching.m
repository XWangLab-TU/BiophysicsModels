function [obj] = Branching(obj,Las17,Mesh,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModActin'));
ip.addRequired('Las17', @(x) isa(x,'ModFreeParticle'));
ip.addRequired('Mesh', @(x) isstruct(x));
ip.addParameter('dt', [], @isnumeric); 
ip.parse(obj,Las17,Mesh,varargin{:});
%--------------------------------------------------------------------------------------------------------
%================================================================================
dt=ip.Results.dt;
if isempty(dt)
    dt=obj.pm.dt;
end
[~,id] = ComMath.getMeshBasedDmin(obj.var.coord,Las17.var.coord,Mesh,'unq1or2',1);
NsuForBr=size(id,1);
angTry=(0:0.1*pi:1.9*pi);
angFix=ones(size(angTry))*obj.pm.AngBr;
nTry=numel(angTry);
coordBest=ones(nTry,3);
for IsuForBr=1:NsuForBr
    iSU=id(IsuForBr,1);
    iBr=obj.var.iSUinBr(iSU);
    rdnNum=rand(1,1);
    if rdnNum<obj.pm.kbr*dt
        for iTry=1:nTry
        coord_new=ComMath.rotAxis([0;0;obj.pm.Dsu],angFix(iTry),angTry(iTry));
        coord_new=ComMath.rotAxis(coord_new,obj.var.angBr(iBr,1),obj.var.angBr(iBr,2));
        coordBest(iTry,:)=coord_new'+obj.var.coord(iSU,:);
        end
        D=sum((Las17.var.coord(id(IsuForBr,2),:)-coordBest).^2,2);
        [~,idBest]=min(D);
        canBr=true;
    else
        canBr=false;
    end
    if canBr==true
        obj.var.coord=[obj.var.coord; coordBest(idBest,:)];
        obj.var.n_coord=obj.var.n_coord+1;
        coord_new=coordBest(idBest,:)-obj.var.coord(iSU,:);
        rXY=sqrt(coord_new(1)^2+coord_new(2)^2);
        if coord_new(1)>0
            obj.var.angBr=[obj.var.angBr;[acos(coord_new(3)/obj.pm.Dsu),pi-acos(coord_new(2)/rXY)]];
        else
            obj.var.angBr=[obj.var.angBr;[acos(coord_new(3)/obj.pm.Dsu),pi+acos(coord_new(2)/rXY)]];
        end
        obj.var.Nbr=obj.var.Nbr+1;
        obj.var.lBr=[obj.var.lBr;1];
        obj.var.capped=[obj.var.capped;false];
        obj.var.iSUinBr=[obj.var.iSUinBr;obj.var.Nbr];
    end
end
%%
% f=figure;
% ComPlot.Sphere(obj.var.coord(1,1),obj.var.coord(1,2),obj.var.coord(1,3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[0 1 0]);
% for iSU=2:4
%     ComPlot.Sphere(obj.var.coord(iSU,1),obj.var.coord(iSU,2),obj.var.coord(iSU,3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[1 0 0]);
% end
% ComPlot.Sphere(0,0,0,50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[0 1 0]);
% ComPlot.Sphere(coord_new(1),coord_new(2),coord_new(3),50,'re_size',obj.pm.Dsu*0.5,'f', f,'col',[0 0 1]);
% xlabel('x');ylabel('y');zlabel('z');
%================================================================================
end