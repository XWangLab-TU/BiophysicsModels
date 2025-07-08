function [m] = ModMemAdapter_ModMembrane8Helfrich(f,m,adp,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('m', @(x) isa(x,'ModMembrane'));
ip.addRequired('adp', @(x) isa(x,'ModMemAdapter'));
ip.addParameter('Cadp', [], @isnumeric);
ip.addParameter('idx', [], @isnumeric);
ip.addParameter('init', [], @islogical);
ip.parse(f,m,adp,varargin{:});
%----------------------------------------------------------------------------------------
idx = ip.Results.idx;
if isempty(idx)
    idx = (1:m.var.n_coord);
end
if ip.Results.init==true
   m.var.f.kH=zeros(m.var.n_coord,1);
   m.var.f.A=zeros(m.var.n_coord,1);
   m.var.f.H=0;
   m.var.f.u_K=zeros(m.var.n_coord,3);
   m.var.f.K=zeros(m.var.n_coord,3);
   m.var.f.dH=zeros(m.var.n_coord,3);
end
Cadp = ip.Results.Cadp;
if isempty(Cadp)
    Cadp=zeros(m.var.n_coord,1);
end
%----------------------------------------------------------------------------------------
%========================================================================================
%%
n_idx = max(size(idx));
for i_ver = 1:n_idx
    i = idx(i_ver);
    c = (m.var.coord(i,:)-m.var.coord(m.var.j_T(i,m.var.j_idx{i}(1:end-1))',:));
    b = m.var.coord(m.var.j_T(i,m.var.j_idx{i}(2:end))',:)-m.var.coord(m.var.j_T(i,m.var.j_idx{i}(1:end-1))',:);
    c_d_c = sum(c.*c,2);
    b_d_b = sum(b.*b,2);
    b_d_c = sum(b.*c,2);
    b1_d_c2 = sum(b.*c(m.var.j_idx{i}(2:end),:),2);
    c2_d_c2 = sum(c(m.var.j_idx{i}(2:end),:).*c(m.var.j_idx{i}(2:end),:),2);
    c1_d_c2 = sum(c.*c(m.var.j_idx{i}(2:end),:),2);
    %-------------------
    A_T_org = 0.5*sqrt(c_d_c.*b_d_b-(b_d_c).^2);
    cosine_abg = [b_d_c./(sqrt(b_d_b).*sqrt(c_d_c)),...
                        -b1_d_c2./(sqrt(b_d_b).*sqrt(c2_d_c2)),...
                        c1_d_c2./(sqrt(c_d_c).*sqrt( c2_d_c2 ) ) ];
    id_tem_all = HeavisideStepContFunc(-cosine_abg,100);
    A_T = 0.5*A_T_org.*id_tem_all(:,3);
    A_T = A_T+0.25*A_T_org.*id_tem_all(:,1);
    A_T = A_T+0.25*A_T_org.*id_tem_all(:,2);
    %-------------------
    Dela = sqrt(c_d_c.*b_d_b-b_d_c.^2);
    Delb = sqrt(c_d_c(m.var.j_idx{i}(2:end)).*b_d_b-b1_d_c2.^2);
    cota = b_d_c./Dela;
    cotb = -b1_d_c2./Delb;
    A_tem = 0.125*(cota+cotb(m.var.j_idx{i}(2:end))).*c_d_c(m.var.j_idx{i}(2:end));
    cos_min = min(cosine_abg,[],2);
    id_tem = HeavisideStepContFunc(cos_min,100);
    A_tem = A_tem.*id_tem;
    A_tem = A_tem+A_T;
    m.var.f.A(i) = sum(A_tem);
    m.var.f.K(i,:) = 0.5/m.var.f.A(i)*sum((cota+cotb(m.var.j_idx{i}(2:end))).*c(m.var.j_idx{i}(2:end),:));
    m.var.f.kH(i) = 0.5*sqrt(sum(m.var.f.K(i,:).*m.var.f.K(i,:),2));
end
%%
%========================================================================================
m.var.f.H = sum((m.var.f.kH-Cadp).^2.*m.var.f.A);
%----------------------------------------------------------------------------------------
end
%========================================================================================
%========================================================================================
function [f] = HeavisideStepContFunc(x,k,varargin)
f = 0.5+0.5.*tanh(k.*x);
end