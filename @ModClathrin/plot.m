function [f] = plot(obj,varargin)
%--------------------------------------------------------------------------
        % plot performs the visulization of @ModClathrin
        % input: 
        % obj - a @ModClathrin object
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('obj', @(x) isa(x,'ModClathrin'));
ip.addParameter('f', [], @isobject);
ip.addParameter('col', [], @isnumeric);
ip.addParameter('FaceAlpha', 1, @isnumeric);
ip.addParameter('LineStyle', 'none', @ischar);
ip.addParameter('LineWidth', 1, @isnumeric);
ip.addParameter('proxy_only', false, @islogical);
ip.addParameter('simple', false, @islogical);
ip.addParameter('iC_iLeg', [], @isnumeric);
ip.parse(obj, varargin{:}); 
%----------------------------------------------------------------------------------------
col=ip.Results.col;
proxy_only=ip.Results.proxy_only;
FaceAlpha=ip.Results.FaceAlpha;
if isempty(col)
    col=rand(obj.var.n_coord,3);
end
if isempty(ip.Results.f)
    f=figure;hold on;
else
    figure(ip.Results.f); hold on;
end
%----------------------------------------------------------------------------------------
%%
if ip.Results.simple==false
[x,y,z] = sphere(10);
for i=1:obj.var.n_coord
            hold on;
            if proxy_only==false
                coord_tem=obj.get_r_from_a(obj.var.coord_org(:,1:3),52,obj.var.O(:,:,i));
            else
                coord_tem=obj.get_r_from_a(obj.var.coord_org(1:24,1:3),24,obj.var.O(:,:,i));
            end
            X=obj.var.coord(i,1)+coord_tem(:,1);
            Y=obj.var.coord(i,2)+coord_tem(:,2);
            Z=obj.var.coord(i,3)+coord_tem(:,3);
            n = max(size(X));
            for i_plot = 1:n
                re_size=obj.var.coord_org(i_plot,4);
                s1=surf(x*re_size+X(i_plot),y*re_size+Y(i_plot),z*re_size+Z(i_plot));
                %set(s1,'facecolor',col_tem,'EdgeColor','none');
%                 set(s1,'facecolor',col(i,:),'EdgeColor','none','facealpha',1);
                set(s1,'facecolor',col(i,:),'FaceLighting','gouraud','EdgeColor','none','facealpha',FaceAlpha); 
            end
            light('Position',[-0.2 -0.2 1],'Style','infinite');
end
else
    iC_iLeg=ip.Results.iC_iLeg;
    for i=1:obj.var.n_coord
            hold on;
            if proxy_only==false
                id_tem=[52,8,16,24,49,50,51];
                coord_tem=obj.get_r_from_a(obj.var.coord_org(id_tem,1:3),7,obj.var.O(:,:,i));
            else
                id_tem=[52,8,16,24];
                coord_tem=obj.get_r_from_a(obj.var.coord_org(id_tem,1:3),4,obj.var.O(:,:,i));
            end
            X=obj.var.coord(i,1)+coord_tem(:,1);
            Y=obj.var.coord(i,2)+coord_tem(:,2);
            Z=obj.var.coord(i,3)+coord_tem(:,3);
            if proxy_only==false
                plot3([X((1)),X((2))],[Y((1)),Y((2))],[Z((1)),Z((2))],'linewidth',2,'color',col(i,:));hold on;
                plot3([X((1)),X((3))],[Y((1)),Y((3))],[Z((1)),Z((3))],'linewidth',2,'color',col(i,:));hold on;
                plot3([X((1)),X((4))],[Y((1)),Y((4))],[Z((1)),Z((4))],'linewidth',2,'color',col(i,:));hold on;
                plot3([X((2)),X((5))],[Y((2)),Y((5))],[Z((2)),Z((5))],'linewidth',0.5,'color',col(i,:));hold on;
                plot3([X((3)),X((6))],[Y((3)),Y((6))],[Z((3)),Z((6))],'linewidth',0.5,'color',col(i,:));hold on;
                plot3([X((4)),X((7))],[Y((4)),Y((7))],[Z((4)),Z((7))],'linewidth',0.5,'color',col(i,:));hold on;
            else
                plot3([X((1)),X((2))],[Y((1)),Y((2))],[Z((1)),Z((2))],'linewidth',2,'color',col(i,:));hold on;
                plot3([X((1)),X((3))],[Y((1)),Y((3))],[Z((1)),Z((3))],'linewidth',2,'color',col(i,:));hold on;
                plot3([X((1)),X((4))],[Y((1)),Y((4))],[Z((1)),Z((4))],'linewidth',2,'color',col(i,:));hold on;
            end
            if ~isempty(iC_iLeg)
                if ismember(i,iC_iLeg(:,1))
                    id_c=iC_iLeg(:,1)==i;    
                    id_leg=iC_iLeg(id_c,2);
                    n_leg=numel(id_leg);
                    for i_leg=1:n_leg
                    X_save=[X((1)),X((id_leg(i_leg)+1))];
                    Y_save=[Y((1)),Y((id_leg(i_leg)+1))];
                    Z_save=[Z((1)),Z((id_leg(i_leg)+1))];
%                     col_tem=rand(1,3);
                    col_tem=[0 1 0];
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    if proxy_only==false
                    X_save=[X(id_leg(i_leg)+1),X((id_leg(i_leg)+4))];
                    Y_save=[Y(id_leg(i_leg)+1),Y((id_leg(i_leg)+4))];
                    Z_save=[Z(id_leg(i_leg)+1),Z((id_leg(i_leg)+4))];
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    plot3([X_save(1),X_save(2)],[Y_save(1),Y_save(2)],[Z_save(1),Z_save(2)],'--','linewidth',4,'color',col_tem);hold on;
                    end
                    end
                end
            end
   end
end
            xlabel('x');ylabel('y');zlabel('z');
            x_tem=xlim;y_tem=ylim;z_tem=zlim;
            lim_tem=[x_tem;y_tem;z_tem];
            lim_tem=[min(lim_tem(:,1)),max(lim_tem(:,2))];
            xlim(lim_tem);ylim(lim_tem);zlim(lim_tem);
%axis equal

