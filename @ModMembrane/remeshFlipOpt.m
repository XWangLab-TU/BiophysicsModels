function [m,remeshed] = remeshFlipOpt(m,j,i_edg,varargin)
%--------------------------------------------------------------------------
        % remeshFlipOpt performs the flipping on @ModMembrane,
        % input: 
        % m - @ModMembrane input
        % j - flipping target [j(1) j(2)]
        % i_edg - index of original edge
        % output:
        % remeshed - true: remesh done, false: remesh gave up
        % optional:
        % see variable arguments
        %   See also remeshSplitOpt,remeshMergeOpt
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = true;
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('j', @(x) isnumeric(x));
ip.addRequired('i_edg', @(x) isnumeric(x));
ip.parse(m,j,i_edg,varargin{:});
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%%
%----------------------------------------------------------------------
remeshed=true;
var_new=m.var;
j = flip.j{i_edg}; i_e = m.var.edge_all(i_edg,:);
%----------------------------------------------------------------------
var_new.edge_all(i_edg,:)=[j(1) j(2)];
%----------------------------------------------------------------------
FaceToAdd=zeros(2,3);
idAll=(1:size(m.var.face_unq,1))';
idTem=sum((m.var.face_unq==i_e(1))+(m.var.face_unq==j(1)),2);
idTem=idAll(idTem==2);
if ismember(i_e(2),m.var.face_unq(idTem(1),:))
    idTem(1,:)=[];
else
    idTem(2,:)=[];
end

nZeroMax=0;
faceAttempt=[i_e(1) j(1) 0; 0 i_e(1) j(1); j(1) 0 i_e(1)];
%                 faceAttempt=[j(1) i_e(1) 0; 0 j(1) i_e(1); i_e(1) 0 j(1)];
faceOrg=m.var.face_unq(idTem,:);
for iA=1:3
    Dface=faceOrg-faceAttempt(iA,:);
    nZero=numel(Dface(Dface==0));
    if nZeroMax<nZero
        nZeroMax=nZero;
    end
end
if nZeroMax==2
    FaceToAdd(1,:)=[i_e(1) j(2) j(1)];
else
    FaceToAdd(1,:)=[i_e(1) j(1) j(2)];
end

idTem=sum((m.var.face_unq==i_e(2))+(m.var.face_unq==j(1)),2);
idTem=idAll(idTem==2);
if ismember(i_e(1),m.var.face_unq(idTem(1),:))
    idTem(1,:)=[];
else
    idTem(2,:)=[];
end

nZeroMax=0;
faceAttempt=[i_e(2) j(1) 0; 0 i_e(2) j(1); j(1) 0 i_e(2)];
faceOrg=m.var.face_unq(idTem,:);
for iA=1:3
    Dface=faceOrg-faceAttempt(iA,:);
    nZero=numel(Dface(Dface==0));
    if nZeroMax<nZero
        nZeroMax=nZero;
    end
end
if nZeroMax==2
    FaceToAdd(2,:)=[i_e(2) j(2) j(1)];
else
    FaceToAdd(2,:)=[i_e(2) j(1) j(2)];
end

idFaceToRem=sum((m.var.face_unq==i_e(1))+(m.var.face_unq==i_e(2)),2);
idFaceToRem=idAll(idFaceToRem==2);
%----------------------------------------------------------------------
var_new.face_unq(idFaceToRem,:)=FaceToAdd;
%----------------------------------------------------------------------
[m.var,topologicalDefect] = m.remeshAddVertex(m.pm,m.var,var_new.coord,var_new.edge_all,var_new.face_unq);
if topologicalDefect==true
    m.failInfo='topologicalDefect';
else
    [m] = getUface(m);
end
%----------------------------------------------------------
%--------------------------------------------------------------------------    
end
