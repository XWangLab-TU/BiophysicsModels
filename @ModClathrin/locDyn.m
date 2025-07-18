function [mod] = locDyn(c,mod,name,idClathrin,idOther,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('name', @(x) iscell(x));
ip.addRequired('idClathrin', @(x) isnumeric(x));
ip.addRequired('idOther', @(x) isnumeric(x));
ip.addParameter('sMax', 0.1, @isnumeric);
ip.addParameter('sMaxCutoff', 1e-6, @isnumeric);
ip.addParameter('dt', 1e-11, @isnumeric);
ip.addParameter('nt', 1000, @isnumeric);
ip.addParameter('update', true, @islogical);
ip.addParameter('addOtherInfo', false, @islogical);
ip.addParameter('initOtherInfo', false, @islogical);
ip.parse(c,mod,name,idClathrin,idOther,varargin{:});
%--------------------------------------------------------------------------
nF=numel(name);
%--------------------------------------------------------------------------
sMax=ip.Results.sMax;
nt=ip.Results.nt;
dt=ip.Results.dt;
sMaxCutoff=ip.Results.sMaxCutoff;
update=ip.Results.update;
addOtherInfo=ip.Results.addOtherInfo;
initOtherInfo=ip.Results.initOtherInfo;
%--------------------------------------------------------------------------
%%
f=mod.TypForce;
%%
VtotPre=inf;
cPre=mod.mod{mod.i_mod.ModClathrin};
for it=1:nt
    if it>1
        addOtherInfo=false;
        initOtherInfo=false;
    end
    ftot=zeros(mod.mod{mod.i_mod.ModClathrin}.var.n_coord,6);
    Vtot=0;
    for iF=1:nF
    [f,Vtem] = f.(name{iF})(mod,'idSub1',idClathrin,'idSub2',idOther,'initOtherInfo',initOtherInfo,'addOtherInfo', addOtherInfo,'update', update);
    Vtot=Vtot+Vtem;
    iC=f.otherInfo.(name{iF}).identifier(mod.i_mod.ModClathrin);
    ftot(idClathrin,:)=ftot(idClathrin,:)+f.int_comp.(name{iF}){iC}(idClathrin,:);
    end
    [mod.mod{mod.i_mod.ModClathrin}] = dynamics.rotation(mod.mod{mod.i_mod.ModClathrin},ftot,dt,0,'idSub',idClathrin,'sMax',sMax);
    if (VtotPre<Vtot)
        mod.mod{mod.i_mod.ModClathrin}=cPre;
        sMax=sMax*0.1;
    else
        cPre=mod.mod{mod.i_mod.ModClathrin};
        VtotPre=Vtot;
    end
    if sMax<sMaxCutoff
        break;
    end
    fprintf('%d   %f\n',it,Vtot);
end
mod.TypForce=f;
end
