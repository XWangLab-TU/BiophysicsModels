function [f] = storedForce(f,Vname,Vpm,Mname,rPara,varargin)
%--------------------------------------------------------------------------
        % storedForce performs the computation of the constant force (1D)
        % resulting from a linearly interpolated potential (1D)
        % input: 
        % f - @TypForce
        % Vname - potential name, i.e. 'VLennardJones'
        % Vpm - according potential parameter
        % Mname - module name, i.e. 'ModMembrane'
        % rPara - [r1,dr,r2], bounds and step
        % output:
        % f - constant force data structure stored in f.int_const.(Mname)
        % see variable arguments
        %   See also randForce
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/17
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('Vname', @(x) ischar(x));
ip.addRequired('Vpm', @(x) isnumeric(x));
ip.addRequired('Mname', @(x) ischar(x));
ip.addRequired('rPara', @(x) isnumeric(x));
ip.addParameter('rsq_std', 0.5, @isnumeric);
ip.addParameter('std_std', [], @isnumeric);
ip.addParameter('std_fold_dr', 10, @isnumeric);
ip.addParameter('plot_or_not', false, @islogical);
ip.addParameter('avoidInfNan', true, @islogical);
ip.parse(f,Vname,Vpm,Mname,rPara,varargin{:});
%--------------------------------------------------------------------------------------------------------
r1=rPara(1);
dr=rPara(2);
r2=rPara(3);
std_fold_dr=ip.Results.std_fold_dr;
std_std=ip.Results.std_std;
if isempty(std_std)
std_std=std_fold_dr*dr;
end
plot_or_not = ip.Results.plot_or_not;
rsq_std = ip.Results.rsq_std;
avoidInfNan=ip.Results.avoidInfNan;
%%
fc =  struct(  ...
               'fn', [],'Vn', [], 'rn', [], 'nr', [], 'in', [],'rg', [],'ig', [],'ng', [],'R', []...
);%fc-force_const; fn-fine scaled force; rn-fine scale; nr-numel of rn; rg-grouped scale; ng-numel of rg; R-stored 'stable' random number, used in ModMembrane
%%
fc.R = random('Stable',1,1,1,0,[10000,1]);
%--------------------------------------------------------------------------------------------------------
%%
fc.rn = (r1+dr:dr:r2-dr);fc.rn = fc.rn';
fc.nr = max(size(fc.rn,1));
Vn=potential.(Vname)([0 0 0],[fc.rn,zeros(fc.nr,2)],Vpm);
%%
i_sad = 1; %saddle point
eps = 0.0;
for i = 2:fc.nr-1
    if ((Vn(i-1) < Vn(i)-eps) && (Vn(i+1) < Vn(i)-eps)) || ((Vn(i-1) > Vn(i)+eps) && (Vn(i+1) > Vn(i)+eps))
        i_sad = [i_sad;i];
    end
end
i_sad = [i_sad;fc.nr];
if plot_or_not
    figure;
plot(fc.rn,Vn); hold on;
plot(fc.rn(i_sad),Vn(i_sad),'*'); hold on;
%ylim([0 8])
end
%%
fc.ig=i_sad;
seg=[fc.ig(1:end-1),fc.ig(2:end)];
n_seg=size(seg,1);
ig_new=fc.ig;
settle=false(n_seg,1);
while numel(settle(settle==false))>0
    min_rsq=inf; std_max=0;
    settle=false(n_seg,1);
for i_seg=1:n_seg
    if seg(i_seg,2)-seg(i_seg,1)>1
    y = Vn(seg(i_seg,1):seg(i_seg,2));
    x = fc.rn(seg(i_seg,1):seg(i_seg,2));
    p = polyfit(x,y,1);
    yfit =  p(1) * x + p(2);
    yresid  = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    std_y=std(y);
    if rsq<min_rsq
        min_rsq=rsq;
    end 
    if std_y>std_max
        std_max=std_y;
    end
    if ((rsq < rsq_std) || (std_y > std_std))
        ig_new=[ig_new;floor(0.5*(seg(i_seg,1)+seg(i_seg,2))+0.5)];
        ig_new=sort(ig_new);
    else
        settle(i_seg) = true;
    end
    else
        settle(i_seg) = true;
    end
end
fc.ig=ig_new;
seg=[fc.ig(1:end-1),fc.ig(2:end)];
n_seg=size(seg,1);
end
%fprintf('membrane constant force settled: %d, %f, %f\n',n_seg, min_rsq, std_max);
%%
fc.ng=size(fc.ig,1);
k=zeros(fc.ng-1,1);
d=zeros(fc.ng-1,1);
Vtem=zeros(fc.ng-1,1);
fc.Vn = zeros(fc.nr,1);
for i = 1:fc.ng-1
y = [Vn(fc.ig(i)),Vn(fc.ig(i+1))];
x = [fc.rn(fc.ig(i)),fc.rn(fc.ig(i+1))];
k(i) = (y(2)-y(1))/(x(2)-x(1)); 
d(i)=(y(2)+y(1) -k(i)*(x(2)+x(1)))*0.5;
Vtem(i)=Vn(fc.ig(i+1));
fc.Vn(fc.ig(i):fc.ig(i+1)-1) = d(i)+k(i)*fc.rn(fc.ig(i):fc.ig(i+1)-1);
end
fc.Vn(end)=Vn(end);

fc.fn = zeros(fc.nr,1);
fc.in = zeros(fc.nr,1);
for i = 2:fc.ng
    fc.fn(fc.ig(i-1)+1:fc.ig(i)) = -k(i-1);
    fc.in(fc.ig(i-1)+1:fc.ig(i)) = i;
end
fc.fn(1)=fc.fn(2);
fc.in(1)=fc.in(2);

fc.rg=fc.rn(fc.ig);

f.int_stored.(Mname)=fc;
if plot_or_not

% plot(fc.rn,Vn); hold on;
% % plot(fc.rn(i_sad),Vn(i_sad),'*');hold on;
% plot(fc.rn(fc.ig),Vn(fc.ig),'o');
for i = 1:fc.ng-1
    plot(fc.rn(fc.ig(i):fc.ig(i+1)),fc.Vn(fc.ig(i):fc.ig(i+1)),'-','linewidth',2); hold on;
end
% subplot(1,2,2);
%plot(fc.rn,f_std); hold on;
%ylim([-70 70]);
end
%==============================================================================
%==============================================================================
end
%==============================================================================

