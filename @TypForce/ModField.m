function [f,Vtot] = ModField(f,M,para,varargin)
%--------------------------------------------------------------------------
        % ModField performs the computation of force fields
        % input: 
        % f - @TypForce
        % M - @Model 
        % para - parameter
        %      i_mod - indicating which module(s) is under field force
        %      type - field type
        %      spec - specific parameter values
        % output:
        % Vtot - total potential in fields
        % optional:
        % see variable arguments
        %   See also TypForce
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2023/11/13
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('M', @(x) isa(x,'model'));
ip.addRequired('para', @(x) isstruct(x));
ip.addParameter('plot_or_not', false, @islogical);
ip.parse(f,M,para,varargin{:});
%----------------------------------------------------------------------------------------
plot_or_not=ip.Results.plot_or_not;
i_mod=para.i_mod;
type=para.type;
spec=para.spec;
n_mod=numel(i_mod);
Vtot=cell(n_mod,1);
typePool={'stripe_z','stripe_z_hole'};
nPool=numel(typePool);
%----------------------------------------------------------------------------------------
%% initiate fields
%----------------------------------------------------------------------------------------
for i=1:n_mod
    f.int_field{i_mod(i)}=zeros(size(M.mod{i_mod(i)}.var.coord));
end
%----------------------------------------------------------------------------------------
%% compute field
%----------------------------------------------------------------------------------------
iPool=0;
for i=1:nPool
    if strcmp(type,typePool{i})
        iPool=i;
        break;
    end
end

for i=1:n_mod
switch iPool
    %----------------------------------------------------------------------
    case 1
        z=M.mod{i_mod(i)}.var.coord(:,3);
        a=spec.steep;
        b=spec.z_lower;
        c=spec.z_upper;
        d=spec.order;
        Vtot{i}=exp(-a*(z-b).^d)+exp(a*(z-c).^d);
        f.int_field{i_mod(i)}=zeros(size(M.mod{i_mod(i)}.var.coord));
        f.int_field{i_mod(i)}(:,3)=a*d*((z-b).^(d-1).*exp(-a*(z-b).^d)-(z-c).^(d-1).*exp(a*(z-c).^d));
        if plot_or_not
            %%
            z=(-10:0.01:10);
            subplot(2,1,1);
            plot(z,exp(-a*(z-b).^d)+exp(a*(z-c).^d)); xlim([min(z), max(z)]); ylim([0 10])
            subplot(2,1,2);
            plot(z, a*d*((z-b).^(d-1).*exp(-a*(z-b).^d)-(z-c).^(d-1).*exp(a*(z-c).^d)) ); xlim([min(z), max(z)]); ylim([-100 100])
        end
    %----------------------------------------------------------------------    
    case 2
        x=M.mod{i_mod(i)}.var.coord(:,1);
        y=M.mod{i_mod(i)}.var.coord(:,2);
        z=M.mod{i_mod(i)}.var.coord(:,3);
        a=spec.steep;
        b=spec.z_lower;
        c=spec.z_upper;
        d=spec.order;
        e=spec.r;
        Vtot{i}=exp(-a*(z-b).^d)+exp(a*(z-c).^d);
        f.int_field{i_mod(i)}=zeros(size(M.mod{i_mod(i)}.var.coord));
        f.int_field{i_mod(i)}(:,3)=a*d*((z-b).^(d-1).*exp(-a*(z-b).^d)-(z-c).^(d-1).*exp(a*(z-c).^d));
        idClear=false(size(M.mod{i_mod(i)}.var.coord));
        idClear(((x.^2+y.^2)<e*e) & (z>0),:)=true;
        f.int_field{i_mod(i)}(idClear)=0;
    %----------------------------------------------------------------------
    otherwise
        error('no supported type found');
end
end
%----------------------------------------------------------------------------------------