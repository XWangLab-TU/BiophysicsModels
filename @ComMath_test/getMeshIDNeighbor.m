function [j] = getMeshIDNeighbor(i,Mesh_range,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('i', @(x) isnumeric(x));
ip.addRequired('Mesh_range', @(x) iscell(x));
ip.parse(i,Mesh_range,varargin{:});
%--------------------------------------------------------------------------
%%
ny=numel(Mesh_range{2});
nz=numel(Mesh_range{3});
%--------------------------------------------------------------------------
%%
j=[i; i+1; i-1; i+nz; i-nz; i+ny*nz; i-ny*nz;... %updown frontback leftright
   i+nz+1; i+nz-1; i-nz+1; i-nz-1;...
   i+ny*nz+1; i+ny*nz-1; i-ny*nz+1; i-ny*nz-1;...
   i+(ny+1)*nz; i+(ny-1)*nz; i-(ny+1)*nz; i-(ny-1)*nz;...
   i+(ny+1)*nz+1; i+(ny+1)*nz-1; i+(ny-1)*nz+1; i+(ny-1)*nz-1; i-(ny+1)*nz+1; i-(ny-1)*nz+1; i-(ny+1)*nz-1; i-(ny-1)*nz-1];
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',[0 0 1]);
%--------------------------------------------------------------------------
%%
% scatter3(mod.Mesh.coord(i,1),mod.Mesh.coord(i,2),mod.Mesh.coord(i,3),'filled','markerfacecolor',[0 0 1]);hold on;
% scatter3(mod.Mesh.coord(i+1,1),mod.Mesh.coord(i+1,2),mod.Mesh.coord(i+1,3),'filled','markerfacecolor',[1 0 1]);hold on;
% scatter3(mod.Mesh.coord(i-1,1),mod.Mesh.coord(i-1,2),mod.Mesh.coord(i-1,3),'filled','markerfacecolor',[1 0 1]);hold on;
% scatter3(mod.Mesh.coord(i+nz,1),mod.Mesh.coord(i+nz,2),mod.Mesh.coord(i+nz,3),'filled','markerfacecolor',[1 0 1]);hold on;
% scatter3(mod.Mesh.coord(i-nz,1),mod.Mesh.coord(i-nz,2),mod.Mesh.coord(i-nz,3),'filled','markerfacecolor',[1 0 1]);hold on;
% j=i+ny*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',[1 0 1]);hold on;
% j=i-ny*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',[1 0 1]);hold on;
% %--------------------------------------------------------------------------
% col=[0 1 1];
% j=i+nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% %--------------------------------------------------------------------------
% col=[1 1 0];
% j=i+ny*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+ny*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-ny*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-ny*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% %--------------------------------------------------------------------------
% col=[0.5 0.5 1];
% j=i+(ny+1)*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+(ny-1)*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny+1)*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny-1)*nz;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% %--------------------------------------------------------------------------
% col=[1 0.5 0.5];
% j=i+(ny+1)*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+(ny+1)*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+(ny-1)*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i+(ny-1)*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny+1)*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny-1)*nz+1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny+1)*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% j=i-(ny-1)*nz-1;
% scatter3(mod.Mesh.coord(j,1),mod.Mesh.coord(j,2),mod.Mesh.coord(j,3),'filled','markerfacecolor',col);hold on;
% %--------------------------------------------------------------------------
% xlim([-5 5]);ylim([-5 5]);zlim([-5 5]);
