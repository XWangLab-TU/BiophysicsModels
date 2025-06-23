function [M] = docking(m1,m2,M, varargin)
%--------------------------------------------------------------------------
        % docking performs the combination of 2 @Modmembrane, this is
        % required for the vesicle fussion example in @BiopysicsApp
        % input: 
        % m1 - a @Modmembrane object
        % m2 - another @Modmembrane object
        % output: 
        % M - the @model object containing the connected membrane
        % optional:
        % see variable arguments
        %    See also getUface
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('m1', @(x) isobject(x));
ip.addRequired('m2', @(x) isobject(x));
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(m1,m2, M,varargin{:});
%--------------------------------------------------------------------------------------------------------
i_mod=M.i_mod.ModMembrane;
%--------------------------------------------------------------------------------------------------------
%%
centr1=mean(m1.var.coord);
centr2=mean(m2.var.coord);
R12=centr2-centr1;
r12=R12/sqrt(sum(R12.^2));
R1=m1.var.coord-centr1;
r1=R1./sqrt(sum(R1.^2,2));
r1_dot_r12=sum(r1.*r12,2);
[~,id_nearest1]=max(r1_dot_r12);
R21=centr1-centr2;
r21=R21/sqrt(sum(R21.^2));
R2=m2.var.coord-centr2;
r2=R2./sqrt(sum(R2.^2,2));
r2_dot_r21=sum(r2.*r21,2);
[~,id_nearest2]=max(r2_dot_r21);
%--------------------------------------------------------------------------------------------------------
n_face1=size(m1.var.face_unq,1);
id_all1=(1:n_face1)';
id_face1=id_all1(sum(m1.var.face_unq-id_nearest1==0,2)==1);
n_face2=size(m2.var.face_unq,1);
id_all2=(1:n_face2)';
id_face2=id_all2(sum(m2.var.face_unq-id_nearest2==0,2)==1);
%--------------------------------------------------------------------------------------------------------
centr_face1=(m1.var.coord(m1.var.face_unq(id_face1,1),:)+m1.var.coord(m1.var.face_unq(id_face1,2),:)+m1.var.coord(m1.var.face_unq(id_face1,3),:))/3;
centr_face2=(m2.var.coord(m2.var.face_unq(id_face2,1),:)+m2.var.coord(m2.var.face_unq(id_face2,2),:)+m2.var.coord(m2.var.face_unq(id_face2,3),:))/3;
n1=size(centr_face1,1);
n2=size(centr_face2,1);
id_face_nearest=[nan,nan];
d_min=inf;
for i1=1:n1
    for i2=1:n2
        d=norm(centr_face1(i1,:)-centr_face2(i2,:));
        if d<d_min
            d_min=d;
            id_face_nearest=[id_face1(i1),id_face2(i2)];
        end
    end
end
%--------------------------------------------------------------------------------------------------------
P = perms(m1.var.face_unq(id_face_nearest(1),:));
X2=m2.var.coord(m2.var.face_unq(id_face_nearest(2),:),1)';
Y2=m2.var.coord(m2.var.face_unq(id_face_nearest(2),:),2)';
Z2=m2.var.coord(m2.var.face_unq(id_face_nearest(2),:),3)';
X1=m1.var.coord(P,1);
X1=reshape(X1,6,3);
Y1=m1.var.coord(P,2);
Y1=reshape(Y1,6,3);
Z1=m1.var.coord(P,3);
Z1=reshape(Z1,6,3);
d_face=sum((X2-X1).^2+(Y2-Y1).^2+(Z2-Z1).^2,2);
[~,id_tem]=min(d_face);
Popt=P(id_tem,:);
%--------------------------------------------------------------------------------------------------------
coord_docked=[0.5*(X2+X1(id_tem,:))',0.5*(Y2+Y1(id_tem,:))',0.5*(Z2+Z1(id_tem,:))'];
m1.var.coord(Popt,:)=coord_docked;
m2.var.coord(m2.var.face_unq(id_face_nearest(2),:),:)=[];
%--------------------------------------------------------------------------------------------------------
m2.var.edge_all=m2.var.edge_all+m1.var.n_coord;
m2.var.face_unq=m2.var.face_unq+m1.var.n_coord;

face_to_replace2=m2.var.face_unq(id_face_nearest(2),:);
id_tem=sum(m2.var.edge_all==m2.var.face_unq(id_face_nearest(2),1),2)...
      +sum(m2.var.edge_all==m2.var.face_unq(id_face_nearest(2),2),2)...
      +sum(m2.var.edge_all==m2.var.face_unq(id_face_nearest(2),3),2);
id_tem=id_tem==2;
m2.var.edge_all(id_tem,:)=[];
m2.var.face_unq(id_face_nearest(2),:)=[];
m1.var.face_unq(id_face_nearest(1),:)=[];

for i=1:3
m2.var.face_unq(m2.var.face_unq==face_to_replace2(i))=Popt(i);
m2.var.edge_all(m2.var.edge_all==face_to_replace2(i))=Popt(i);
end
face_to_replace2=sort(face_to_replace2);
for i=1:3
    id_tem=m2.var.edge_all>=face_to_replace2(i);
    m2.var.edge_all(id_tem)=m2.var.edge_all(id_tem)-1;
    id_tem=m2.var.face_unq>=face_to_replace2(i);
    m2.var.face_unq(id_tem)=m2.var.face_unq(id_tem)-1;
    face_to_replace2=face_to_replace2-1;
end
%--------------------------------------------------------------------------------------------------------
%%
coord_new=[m1.var.coord;m2.var.coord];
edge_all_new=[m1.var.edge_all;m2.var.edge_all];
face_unq_new=[m1.var.face_unq;m2.var.face_unq];
[M.mod{i_mod}.var] = M.mod{i_mod}.remeshAddVertex(M.mod{i_mod}.pm,M.mod{i_mod}.var,coord_new,edge_all_new,face_unq_new);
[M.mod{i_mod}] = getUface(M.mod{i_mod});
%--------------------------------------------------------------------------------------------------------
%%
% fig=figure('units','normalized','outerposition',[0 0 1 1]); 
% plot(m1,'linestyle','-','f',fig,'FaceAlpha',1);
% plot(m2,'linestyle','-','f',fig,'FaceAlpha',1);
%--------------------------------------------------------------------------------------------------------
%%
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
% plot(mod.mod{mod.i_mod.ModMembrane},'linestyle','-','f',fig,'FaceAlpha',1);
end
%==========================================================================
