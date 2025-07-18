function [c,changed] = link(c,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @(x) isa(x,'ModClathrin'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('r_min', [], @isnumeric);
ip.parse(c,varargin{:});
%--------------------------------------------------------------------------
r_min=ip.Results.r_min;
if isempty(r_min)
%  r_min=norm(mod.mod{1}.var.coord_org(1,1:3)-mod.mod{1}.var.coord_org(2,1:3))*10;
   r_min=c.pm.r_min_for_link;
end
%----------------------------------------------------------------------------------------
%mex_avail=ip.Results.mex_avail;
mex=ip.Results.mex;
if isempty(mex)
    mex_avail=false;
else
    mex_avail=true;
end
%----------------------------------------------------------------------------------------
changed=false;
[waist] = getWaist(c);
% [CtrlPt] = getCtrlPt(mod.mod{i_mod});
% u=zeros(mod.mod{i_mod}.var.n_coord,3);
% for i_c=1:mod.mod{i_mod}.var.n_coord
% u(i_c,:)=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c));
% end
%%
id_all=(1:3);
for i_c1=1:c.var.n_coord-1
%     disp(i_c1);
    if ismember(0,c.var.connect(:,1,i_c1))
    for i_c2=i_c1+1:c.var.n_coord
%         disp(i_c2);
        %cos_tem=dot(u(i_c1,:),u(i_c2,:));%norm(u(i_c1,:)) and norm(u(i_c2,:)) = 1
        if (~ismember(i_c2,c.var.connect(:,1,i_c1))) ...
          &&(ismember(0,c.var.connect(:,1,i_c2)))
%--------------------------------------------------------------------------  
            id1=id_all(c.var.connect(:,1,i_c1)==0);
            id2=id_all(c.var.connect(:,1,i_c2)==0);
            d1=sqrt(sum((waist(id1,:,i_c1)-c.var.coord(i_c2,:)).^2,2));
            d2=sqrt(sum((waist(id2,:,i_c2)-c.var.coord(i_c1,:)).^2,2));
            id1=id1(d1<r_min);
            id2=id2(d2<r_min);
            if (~isempty(id1)) && (~isempty(id2))
                %d_Ctrl=sqrt(sum((CtrlPt(id1,:,i_c1)-CtrlPt(id2,:,i_c2)).^2,2));
                %if d_Ctrl<r_min
%                 disp('yeah');pause;
                connect_save=c.var.connect;
                c.var.connect(id1(1),1,i_c1)=i_c2;
                c.var.connect(id1(1),2,i_c1)=id2(1);
                c.var.connect(id2(1),1,i_c2)=i_c1;
                c.var.connect(id2(1),2,i_c2)=id1(1);
                junction=[i_c1 id1(1);i_c2 id2(1)];
%                 disp('before');
%                 disp(junction);
%                 pause;
                [topology,n_ring] = getTopology(c,junction);
%                 disp('after');
                if (topology==true) && (n_ring>4)
                    fprintf('new connection added: %d,%d; %d,%d\n',i_c1,id1(1),i_c2,id2(1));
                    changed=true;
                else
                    c.var.connect=connect_save;
                end
                
                %end
            end

%--------------------------------------------------------------------------
        end
    end
    end
end

