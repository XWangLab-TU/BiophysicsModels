function [dyn] = Preparation(dyn,M,Fname,varargin)
%--------------------------------------------------------------------------
        % Preparation performs the preparation for dynamics to determin
        % which type of motions and match forces to modules
        % input: 
        % dyn - a @dynamics object
        % M - @model object including all of the modules in dynamics
        % Fname - name of forces to be included, e.g. Fname={'ModMembrane','ModMembrane_ModSubstrate'};
        % optional:
        % see variable arguments
        %   See also TimeEval
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dyn', @(x) isa(x,'dynamics'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('Fname', @(x) iscell(x));
ip.addParameter('field_or_not', [], @islogical); 
ip.addParameter('field_para', [], @iscell); %add field forces parameter
%e.g.: para{1}=struct('i_mod',3,'type','stripe_z','spec',struct('z_lower',-1,'z_upper',1,'order',1,'steep',5));
ip.parse(dyn,M,Fname,varargin{:});
%--------------------------------------------------------------------------
nF=numel(Fname);
field_or_not=ip.Results.field_or_not;
field_para=ip.Results.field_para;
%--------------------------------------------------------------------------
dyn.ifVarDt=false(M.n_mod,1);
dyn.needBreakOff=false(M.n_mod,1);
dyn.updateModForce=false(M.n_mod,1);
dyn.needFollow=false(M.n_mod,1);
dyn.needField=false(M.n_mod,1);
for im=1:M.n_mod
    for ip=1:numel(M.mod{im}.prop)
        if strcmp(M.mod{im}.prop{ip},'varDt')
             dyn.ifVarDt(im)=true;
        elseif strcmp(M.mod{im}.prop{ip},'needBreakOff')
            dyn.needBreakOff(im)=true;
        elseif strcmp(M.mod{im}.prop{ip},'forceRelateMod')
            dyn.updateModForce(im)=true;
        elseif strcmp(M.mod{im}.prop{ip},'needFollow')
            if isempty(M.mod{im}.var.follow)
                error('need to assign follow info');
            end
            dyn.needFollow(im)=true;
        end
    end
end
if field_or_not  
    dyn.field=cell(M.n_mod,1);
    n_field=numel(field_para);
    for i_field=1:n_field
        dyn.field{field_para{i_field}.i_mod}=field_para{i_field};
    end
end
%--------------------------------------------------------------------------
%%
dyn.dynTyp=-ones(M.n_mod,1); %0: Langevin only; 1: rotation; -1: fixed
for im=1:M.n_mod
    if strcmp(M.mod{im}.prop{1},'Particle')
        dyn.dynTyp(im)=0;
    elseif strcmp(M.mod{im}.prop{1},'RigidBody')
        dyn.dynTyp(im)=1;
    elseif strcmp(M.mod{im}.prop{1},'Fixed')
        dyn.dynTyp(im)=-1;
    end
end
%first col: 0: single force; >0: col number of F assigned to a mod
%second cols: corresponding number of F
%rows: muliple F to one mod 
dyn.iFon=false(nF,1);
dyn.matchFtoMod=cell(M.n_mod,1); 
for im=1:M.n_mod
    dyn.matchFtoMod{im}=[];
    for iF=1:nF
        if strcmp(Fname{iF},M.name{im})
            dyn.matchFtoMod{im}=[dyn.matchFtoMod{im};[1 iF]];
            dyn.iFon(iF)=true;
        else
            idTem=Fname{iF}=='_';
            SidTem=sum(idTem);
            idAll=1:numel(Fname{iF});
            if SidTem==1
                idBar=idAll(idTem);
                if strcmp(Fname{iF}(1:idBar-1),M.name{im})
                    dyn.matchFtoMod{im}=[dyn.matchFtoMod{im};[1 iF]];
                    dyn.iFon(iF)=true;
                elseif strcmp(Fname{iF}(idBar+1:end),M.name{im})
                    dyn.matchFtoMod{im}=[dyn.matchFtoMod{im};[2 iF]];
                    dyn.iFon(iF)=true;
                end
            elseif SidTem>1
                disp(Fname{iF});
                error('type of interaction currently not supported');
            end                
        end
    end
end

end
