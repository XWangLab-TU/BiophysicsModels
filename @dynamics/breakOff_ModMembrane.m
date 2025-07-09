function [M,breakOffInfo] = breakOff_ModMembrane(dyn,M,varargin)
%--------------------------------------------------------------------------
        % breakOff_ModMembrane performs the breakoff operation for
        % @ModMembrane during dynamics, here the remeshing for @ModMembrane
        % and the lipid transfer under diffusion barrier for tetherPull in
        % the @BiophysicsApp is done
        % input: 
        % dyn - a @dynamics object
        % M - @model object including all of the modules in dynamics
        % output:
        % breakOff - true, remeshing needed; false, remeshing not needed
        % optional:
        % see variable arguments
        %   See also TimeEval, breakOff_ModClathrin
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------  
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dyn', @(x) isa(x,'dynamics'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.parse(dyn,M,varargin{:});
%--------------------------------------------------------------------------
%========================================================================================================================== 
%--------------------------------------------------------------------------
d = (M.mod{M.i_mod.ModMembrane}.var.coord(M.mod{M.i_mod.ModMembrane}.var.edge_all(M.mod{M.i_mod.ModMembrane}.var.id_on_edg,2),:) ...
    -M.mod{M.i_mod.ModMembrane}.var.coord(M.mod{M.i_mod.ModMembrane}.var.edge_all(M.mod{M.i_mod.ModMembrane}.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
%--------------------------------------------------------------------------
breakOffInfo=struct('breakOff',false,'id',[]);
%--------------------------------------------------------------------------
Vpm=M.mod{M.i_mod.ModMembrane}.pm.Vdh;

id_tem1 = r<Vpm.rl_min;  id_tem2 = r>Vpm.rl_max; 
n_tem1 = length(id_tem1(id_tem1));
n_tem2 = length(id_tem2(id_tem2));
if (n_tem1 == 0) && (n_tem2 == 0)
    breakOff=false;
else
    breakOff=true;
    breakOffInfo.id=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
    breakOffInfo.id=[breakOffInfo.id(id_tem1);breakOffInfo.id(id_tem2)];
end
breakOffInfo.breakOff=breakOff;

if breakOff==true
% [M] = ModMembrane(M.TypChemistry,M);
if M.mod{M.i_mod.ModMembrane}.pm.remeshScheme==0
    [M,~] = remesh(M.mod{M.i_mod.ModMembrane},M);
elseif M.mod{M.i_mod.ModMembrane}.pm.remeshScheme==1
    [M,~] = remeshCSbased(M.mod{M.i_mod.ModMembrane},M);
else
    error('remeshScheme not right, 0-physical; 1-CS-based');
end
if isinf(M.mod{M.i_mod.ModMembrane}.pm.kDiff)
%instant diffusion of tension
    M.mod{M.i_mod.ModMembrane}.var.dens=642/M.mod{M.i_mod.ModMembrane}.var.n_coord*ones(M.mod{M.i_mod.ModMembrane}.var.n_coord,1);
else
    %diffusion of tension, define a region with diffusion barrier
    %that does not allow diffusion
    xb=M.mod{M.i_mod.ModMembrane}.pm.xbDiff;
    r_thre=M.mod{M.i_mod.ModMembrane}.pm.rDiff;
    r=sqrt(sum(M.mod{M.i_mod.ModMembrane}.var.coord.^2,2));
    idMid=(abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,1))<xb) & (abs(M.mod{M.i_mod.ModMembrane}.var.coord(:,2)./r)>r_thre);
    idAll=(1:M.mod{M.i_mod.ModMembrane}.var.n_coord);
    idMid=idAll(idMid);
    id_on_edg=(1:M.mod{M.i_mod.ModMembrane}.var.n_edg)';
    id_on_edg(sum(ismember(M.mod{M.i_mod.ModMembrane}.var.edge_all,idMid),2)>0)=[];
    [M.mod{M.i_mod.ModMembrane}] = densTran(M.mod{M.i_mod.ModMembrane},M.mod{M.i_mod.ModMembrane}.pm.kDiff,'nt',50,'id_on_edg',id_on_edg);
end

end