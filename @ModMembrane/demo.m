function [f] = demo(varargin)
%--------------------------------------------------------------------------
        % Author: Xinxin Wang
        % email: wangxinxin8627@gmail.com
        % date: 2025/08/05
%--------------------------------------------------------------------------   
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('plot_or_not', true, @islogical);
ip.parse(varargin{:});
%----------------------------------------------------------------------------------------
f=figure;
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
n=10;
V=zeros(n,2);
r=zeros(n,1);
for l0=1:n
m=ModMembrane(true,3,0,'unit',u,'l0',l0);
r(l0)=max(m.var.coord(:,1));
V(l0,1)=4/3*pi*r(l0)^3;
V(l0,2)=sum(Volume(m));
end
figure(f);
plot(r,V(:,1),'.-'); xlabel('r'); ylabel('V'); title('Volume: threoretical vs computational (Volume func.)'); hold on;
plot(r,V(:,2),'.-'); 
%----------------------------------------------------------------------------------------
end


