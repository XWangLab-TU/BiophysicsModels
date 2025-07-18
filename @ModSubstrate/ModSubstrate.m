classdef ModSubstrate
%==========================================================================
%--------------------------------------------------------------------------
        % ModSubstrate stores the data structure of external control points
        % to mechanically manipulate @ModMembrane, force is defined in
        % @TypForce as spring force
        %   See also ModMembrane, model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
%========================================================================== 
properties
    var %variables
    pm %parameters
    prop %descriptions
    failInfo
end
%==========================================================================
%==========================================================================    
methods
    function obj = ModSubstrate(dr, i_case) 
%--------------------------------------------------------------------------
        % ModSubstrate is the initiation function for @ModSubstrate to setup
        % parameters and variables
        % input: 
        % dr - spatial step length
        % i_case - points to arrays reflecting various substrate cases:
        % 0,1-multi-protrusion; 2-filopodia; 3-lamellipodia; 4-retraction; 5-blebb; 6-flat;
        %   See also TypForce, ModMembrane
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/24
%--------------------------------------------------------------------------
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.addRequired('dr', @(x) isnumeric);
        ip.addRequired('i_case', @(x) isnumeric);
        %             ip.parse(mat_int,dir_mod,varargin{:});
%--------------------------------------------------------------------------
obj.pm.dr=dr;
obj.prop = {'Fixed'};
obj.failInfo='';
%--------------------------------------------------------------------------    
%%
icosphere1 =[...
   -1.0000         0         0;
   -0.8507         0   -0.5257;
   -0.8507         0    0.5257;
   -0.8090   -0.5000   -0.3090;
   -0.8090   -0.5000    0.3090;
   -0.8090    0.5000   -0.3090;
   -0.8090    0.5000    0.3090;
   -0.5257   -0.8507         0;
   -0.5257    0.8507         0;
   -0.5000   -0.3090   -0.8090;
   -0.5000   -0.3090    0.8090;
   -0.5000    0.3090   -0.8090;
   -0.5000    0.3090    0.8090;
   -0.3090   -0.8090   -0.5000;
   -0.3090   -0.8090    0.5000;
   -0.3090    0.8090   -0.5000;
   -0.3090    0.8090    0.5000;
         0   -1.0000         0;
         0   -0.5257   -0.8507;
         0   -0.5257    0.8507;
         0         0   -1.0000;
         0         0    1.0000;
         0    0.5257   -0.8507;
         0    0.5257    0.8507;
         0    1.0000         0;
    0.3090   -0.8090   -0.5000;
    0.3090   -0.8090    0.5000;
    0.3090    0.8090   -0.5000;
    0.3090    0.8090    0.5000;
    0.5000   -0.3090   -0.8090;
    0.5000   -0.3090    0.8090;
    0.5000    0.3090   -0.8090;
    0.5000    0.3090    0.8090;
    0.5257   -0.8507         0;
    0.5257    0.8507         0;
    0.8090   -0.5000   -0.3090;
    0.8090   -0.5000    0.3090;
    0.8090    0.5000   -0.3090;
    0.8090    0.5000    0.3090;
    0.8507         0   -0.5257;
    0.8507         0    0.5257;
    1.0000         0         0];
icosphere2 =[...
   -1.0000         0         0;
   -0.9619         0   -0.2733;
   -0.9619         0    0.2733;
   -0.9511   -0.2629   -0.1625;
   -0.9511   -0.2629    0.1625;
   -0.9511    0.2629   -0.1625;
   -0.9511    0.2629    0.1625;
   -0.8627   -0.2599   -0.4339;
   -0.8627   -0.2599    0.4339;
   -0.8627    0.2599   -0.4339;
   -0.8627    0.2599    0.4339;
   -0.8507   -0.5257         0;
   -0.8507         0   -0.5257;
   -0.8507         0    0.5257;
   -0.8507    0.5257         0;
   -0.8090   -0.5000   -0.3090;
   -0.8090   -0.5000    0.3090;
   -0.8090    0.5000   -0.3090;
   -0.8090    0.5000    0.3090;
   -0.7020   -0.1606   -0.6938;
   -0.7020   -0.1606    0.6938;
   -0.7020    0.1606   -0.6938;
   -0.7020    0.1606    0.6938;
   -0.6938   -0.7020   -0.1606;
   -0.6938   -0.7020    0.1606;
   -0.6938    0.7020   -0.1606;
   -0.6938    0.7020    0.1606;
   -0.6882   -0.4253   -0.5878;
   -0.6882   -0.4253    0.5878;
   -0.6882    0.4253   -0.5878;
   -0.6882    0.4253    0.5878;
   -0.5878   -0.6882   -0.4253;
   -0.5878   -0.6882    0.4253;
   -0.5878    0.6882   -0.4253;
   -0.5878    0.6882    0.4253;
   -0.5257   -0.8507         0;
   -0.5257         0   -0.8507;
   -0.5257         0    0.8507;
   -0.5257    0.8507         0;
   -0.5000   -0.3090   -0.8090;
   -0.5000   -0.3090    0.8090;
   -0.5000    0.3090   -0.8090;
   -0.5000    0.3090    0.8090;
   -0.4339   -0.8627   -0.2599;
   -0.4339   -0.8627    0.2599;
   -0.4339    0.8627   -0.2599;
   -0.4339    0.8627    0.2599;
   -0.4253   -0.5878   -0.6882;
   -0.4253   -0.5878    0.6882;
   -0.4253    0.5878   -0.6882;
   -0.4253    0.5878    0.6882;
   -0.3090   -0.8090   -0.5000;
   -0.3090   -0.8090    0.5000;
   -0.3090    0.8090   -0.5000;
   -0.3090    0.8090    0.5000;
   -0.2733   -0.9619         0;
   -0.2733    0.9619         0;
   -0.2629   -0.1625   -0.9511;
   -0.2629   -0.1625    0.9511;
   -0.2629    0.1625   -0.9511;
   -0.2629    0.1625    0.9511;
   -0.2599   -0.4339   -0.8627;
   -0.2599   -0.4339    0.8627;
   -0.2599    0.4339   -0.8627;
   -0.2599    0.4339    0.8627;
   -0.1625   -0.9511   -0.2629;
   -0.1625   -0.9511    0.2629;
   -0.1625    0.9511   -0.2629;
   -0.1625    0.9511    0.2629;
   -0.1606   -0.6938   -0.7020;
   -0.1606   -0.6938    0.7020;
   -0.1606    0.6938   -0.7020;
   -0.1606    0.6938    0.7020;
         0   -1.0000         0;
         0   -0.8507   -0.5257;
         0   -0.8507    0.5257;
         0   -0.5257   -0.8507;
         0   -0.5257    0.8507;
         0   -0.2733   -0.9619;
         0   -0.2733    0.9619;
         0         0   -1.0000;
         0         0    1.0000;
         0    0.2733   -0.9619;
         0    0.2733    0.9619;
         0    0.5257   -0.8507;
         0    0.5257    0.8507;
         0    0.8507   -0.5257;
         0    0.8507    0.5257;
         0    1.0000         0;
    0.1606   -0.6938   -0.7020;
    0.1606   -0.6938    0.7020;
    0.1606    0.6938   -0.7020;
    0.1606    0.6938    0.7020;
    0.1625   -0.9511   -0.2629;
    0.1625   -0.9511    0.2629;
    0.1625    0.9511   -0.2629;
    0.1625    0.9511    0.2629;
    0.2599   -0.4339   -0.8627;
    0.2599   -0.4339    0.8627;
    0.2599    0.4339   -0.8627;
    0.2599    0.4339    0.8627;
    0.2629   -0.1625   -0.9511;
    0.2629   -0.1625    0.9511;
    0.2629    0.1625   -0.9511;
    0.2629    0.1625    0.9511;
    0.2733   -0.9619         0;
    0.2733    0.9619         0;
    0.3090   -0.8090   -0.5000;
    0.3090   -0.8090    0.5000;
    0.3090    0.8090   -0.5000;
    0.3090    0.8090    0.5000;
    0.4253   -0.5878   -0.6882;
    0.4253   -0.5878    0.6882;
    0.4253    0.5878   -0.6882;
    0.4253    0.5878    0.6882;
    0.4339   -0.8627   -0.2599;
    0.4339   -0.8627    0.2599;
    0.4339    0.8627   -0.2599;
    0.4339    0.8627    0.2599;
    0.5000   -0.3090   -0.8090;
    0.5000   -0.3090    0.8090;
    0.5000    0.3090   -0.8090;
    0.5000    0.3090    0.8090;
    0.5257   -0.8507         0;
    0.5257         0   -0.8507;
    0.5257         0    0.8507;
    0.5257    0.8507         0;
    0.5878   -0.6882   -0.4253;
    0.5878   -0.6882    0.4253;
    0.5878    0.6882   -0.4253;
    0.5878    0.6882    0.4253;
    0.6882   -0.4253   -0.5878;
    0.6882   -0.4253    0.5878;
    0.6882    0.4253   -0.5878;
    0.6882    0.4253    0.5878;
    0.6938   -0.7020   -0.1606;
    0.6938   -0.7020    0.1606;
    0.6938    0.7020   -0.1606;
    0.6938    0.7020    0.1606;
    0.7020   -0.1606   -0.6938;
    0.7020   -0.1606    0.6938;
    0.7020    0.1606   -0.6938;
    0.7020    0.1606    0.6938;
    0.8090   -0.5000   -0.3090;
    0.8090   -0.5000    0.3090;
    0.8090    0.5000   -0.3090;
    0.8090    0.5000    0.3090;
    0.8507   -0.5257         0;
    0.8507         0   -0.5257;
    0.8507         0    0.5257;
    0.8507    0.5257         0;
    0.8627   -0.2599   -0.4339;
    0.8627   -0.2599    0.4339;
    0.8627    0.2599   -0.4339;
    0.8627    0.2599    0.4339;
    0.9511   -0.2629   -0.1625;
    0.9511   -0.2629    0.1625;
    0.9511    0.2629   -0.1625;
    0.9511    0.2629    0.1625;
    0.9619         0   -0.2733;
    0.9619         0    0.2733;
    1.0000         0         0];
%%
%--------------------------------------------------------------------------
        switch i_case
            case 0
                r_tem=4.5;
                ver_tem = icosphere2;
                obj.var.coord = ver_tem*r_tem;
                id_tem = [1 162 89 74 82 81];
                obj.var.coord(id_tem,:) = ver_tem(id_tem,:).*(r_tem+2);
%--------------------------------------------------------------------------
            case 1
                ver_tem = icosphere2;
                r_tem = 6.5;               
                obj.var.coord = ver_tem*r_tem;
                id_tem = [1 162 89 74 82 81];
                obj.var.coord(id_tem,:) = ver_tem(id_tem,:).*(r_tem-2);
%--------------------------------------------------------------------------
            case 2
                ver_tem = icosphere1;
                r_tem = 6.5;
                obj.var.coord = [];
                obj.var.coord(:,1) = ver_tem(:,1)*r_tem;
                obj.var.coord(:,2) = ver_tem(:,2)*r_tem;
                obj.var.coord(:,3) = ver_tem(:,3)*r_tem;
                d_tem=10;
                [~,id_tem] = min(vecnorm(obj.var.coord-[r_tem,0,0],2,2));
                obj.var.coord(id_tem,1) = obj.var.coord(id_tem,1)+d_tem;
                [~,id_tem] = min(vecnorm(obj.var.coord-[r_tem*0.5*sqrt(2),r_tem*0.5*sqrt(2),0],2,2));
                obj.var.coord(id_tem,1) = obj.var.coord(id_tem,1)+d_tem*0.5*sqrt(2);
                obj.var.coord(id_tem,2) = obj.var.coord(id_tem,2)+d_tem*0.5*sqrt(2);
                [~,id_tem] = min(vecnorm(obj.var.coord-[-r_tem*0.5*sqrt(2),r_tem*0.5*sqrt(2),0],2,2));
                obj.var.coord(id_tem,1) = obj.var.coord(id_tem,1)-d_tem*0.5*sqrt(2);
                obj.var.coord(id_tem,2) = obj.var.coord(id_tem,2)+d_tem*0.5*sqrt(2);
%--------------------------------------------------------------------------
            case 3
                d_tem = 10;
                phi = (0.1:0.1:pi-0.1)'+1.5*pi;
                ver_tem = [ d_tem*cos(phi), d_tem*sin(phi), zeros(max(size(phi)),1)];
                ver_tem2 = icosphere1;
                r_tem = 6.5;
                ver_tem2 = ver_tem2*r_tem;
                id_tem = (abs(ver_tem2(:,3))< 0.5) & (ver_tem2(:,1) > 0);
                ver_tem2(id_tem,:) = [];
                obj.var.coord = [ver_tem;ver_tem2];
%--------------------------------------------------------------------------
            case 4                
                ver_tem = icosphere1;
                ver_tem1 = icosphere1;
                r_tem1 = 3;
                d_tem = 0;
                ver_tem1 = ver_tem1*r_tem1;
                id_tem = ver_tem1(:,1) < d_tem;
                ver_tem1(id_tem,:) = [];
                r_tem = 6.5;
                ver_tem = ver_tem*r_tem;
                id_tem = ver_tem(:,1) < -5.5;
                ver_tem(id_tem,:) = [];
                ver_tem1(:,1) = ver_tem1(:,1)-5.5+d_tem;
                ver_tem = [ver_tem;ver_tem1];
                obj.var.coord = ver_tem;
%--------------------------------------------------------------------------
            case 5
                ver_tem = icosphere2;
                ver_tem1 = icosphere1;
                r_tem1 = 3;
                d_tem = -1;
                ver_tem1 = ver_tem1*r_tem1;
                id_tem = ver_tem1(:,1)< d_tem;
                ver_tem1(id_tem,:) = [];
                r_tem = 6.5;
                ver_tem = ver_tem*r_tem;
                id_tem = ver_tem(:,1)> 5.5;
                ver_tem(id_tem,:) = [];
                ver_tem1(:,1) = ver_tem1(:,1)+5.5-d_tem;
                ver_tem = [ver_tem;ver_tem1];
                obj.var.coord = ver_tem;
%--------------------------------------------------------------------------
            case 6
                n_tem = 10;
                ver_tem = zeros(n_tem*n_tem,3);
                ver_tem(:,3) = -n_tem/4;
                for i=1:n_tem
                    ver_tem((i-1)*n_tem+1:i*n_tem,1)=i;
                    ver_tem((i-1)*n_tem+1:i*n_tem,2)=1:n_tem;
                end
                ver_tem(:,1) = ver_tem(:,1) - n_tem/2;
                ver_tem(:,2) = ver_tem(:,2) - n_tem/2;
                %                 ver_tem((n_tem-1)*n_tem+1:i*n_tem,1) = ver_tem((n_tem-1)*n_tem+1:i*n_tem,1)+(rand(n_tem,1)-0.5)*1;
                %                 ver_tem((0)*n_tem+1:1*n_tem,1) = ver_tem((0)*n_tem+1:1*n_tem,1)+(rand(n_tem,1)-0.5)*1;
                %                 ver_tem(n_tem:n_tem:n_tem*n_tem,2) = ver_tem(n_tem:n_tem:n_tem*n_tem,2)+(rand(n_tem,1)-0.5)*1;
                %                 ver_tem(1:n_tem:n_tem*(n_tem-1)+1,2) = ver_tem(1:n_tem:n_tem*(n_tem-1)+1,2)+(rand(n_tem,1)-0.5)*1;
                obj.var.coord = ver_tem*1;
%--------------------------------------------------------------------------
            case 7
                obj.var.coord = [0 0 8];
            case 8
                ver_tem = icosphere1;
                r_tem = 3;               
                obj.var.coord = ver_tem*r_tem;
                id_tem = randsample(size(ver_tem,1), 3);
                obj.var.coord(id_tem,:) = ver_tem(id_tem,:).*(r_tem+1);
            case 9
                d_ang=0.1*pi;
                ang_tem = (d_ang:d_ang:2*pi)';
                r_tem = 7.5;               
                obj.var.coord = [cos(ang_tem)*r_tem, sin(ang_tem)*r_tem, -2*ones(size(ang_tem))];
            %--------------------------------------------------------------------------
            case 10
                ver_tem = icosphere2;
                r_tem = 6.5;
                obj.var.coord = [];
                obj.var.coord(:,1) = ver_tem(:,1)*r_tem;
                obj.var.coord(:,2) = ver_tem(:,2)*r_tem;
                obj.var.coord(:,3) = ver_tem(:,3)*r_tem;
                d_tem=4;
                id_tem = 162;
                obj.var.coord(id_tem,1) = obj.var.coord(id_tem,1)+d_tem;
                [~,id_sort]=sort(obj.var.coord(:,1));
                obj.var.coord(id_sort(end-6:end-1),:)=[];
            %--------------------------------------------------------------------------
            case 11
                ver_tem = icosphere2*3;
                [~,id_tem]=sort(-ver_tem(:,3));
                ver_tem=ver_tem(id_tem(1:40),:)-[0 0 -0.6];
                
                obj.var.coord=ver_tem;
%                 obj.var.coord = [0 0 0];
%                 d_tem=2;
%                 obj.var.coord = [obj.var.coord;[d_tem d_tem 4.]];
%                 obj.var.coord = [obj.var.coord;[d_tem -d_tem 4.]];
%                 obj.var.coord = [obj.var.coord;[-d_tem d_tem 4.]];
%                 obj.var.coord = [obj.var.coord;[-d_tem -d_tem 4.]];
            %--------------------------------------------------------------------------
            case 12
                %%
%                 ang_tem=(0.05*pi:0.05*pi:2*pi)';
%                 r_tem = 6.5;
%                 n_ang=numel(ang_tem);
%                 d=5;
%                 ver_tem = [ones(n_ang,1)*d,sin(ang_tem)*r_tem,cos(ang_tem)*r_tem];
%                 ver_tem = [ver_tem;[-ones(n_ang,1)*d,sin(ang_tem)*r_tem,cos(ang_tem)*r_tem]];
                ver_tem = icosphere2;
                ver_tem(ver_tem(:,1)<-0.01,:)=[];
                d=0.1;
                ver_tem(:,1)=ver_tem(:,1)+d;
                ver_tem=[ver_tem;[-ver_tem(:,1),ver_tem(:,2),ver_tem(:,3)]];            
                r_tem = 6.5;
                ver_tem=ver_tem*r_tem;
                
%                 scatter3(ver_tem(:,1),ver_tem(:,2),ver_tem(:,3));
                obj.var.coord=ver_tem;
            %--------------------------------------------------------------------------
            case 13
                d=6.6;
                obj.var.coord=[-0.5*sqrt(2)*d,0.5*sqrt(2)*d,0;0.5*sqrt(2)*d,0.5*sqrt(2)*d,0;...
                               -0.5*sqrt(2)*d,-0.5*sqrt(2)*d,0;0.5*sqrt(2)*d,-0.5*sqrt(2)*d,0;0,0,d;0,0,-d;];
                obj.var.coord(1:2,:)=obj.var.coord(1:2,:)*3;
%                 obj.var.coord=[d,0,0;-d,0,0;0,d,0;0,-d,0];
%--------------------------------------------------------------------------
            case 14
                r_tem=6.6;
                nG=5;
                theta=linspace(0,2*pi,nG);
                phi=linspace(0,pi,nG);
                [theta,phi]=meshgrid(theta,phi);
                obj.var.coord=[r_tem*reshape(sin(phi).*cos(theta),[nG^2,1]),...
                               r_tem*reshape(sin(phi).*sin(theta),[nG^2,1]),...
                               r_tem*reshape(cos(phi),[nG^2,1])];
                obj.var.coord=[obj.var.coord; r_tem*0.5*sqrt(2),r_tem*0.5*sqrt(2),0;...
                                              -r_tem*0.5*sqrt(2),r_tem*0.5*sqrt(2),0;...
                                              r_tem*0.5*sqrt(2),-r_tem*0.5*sqrt(2),0;...
                                              -r_tem*0.5*sqrt(2),-r_tem*0.5*sqrt(2),0];    
%                 obj.var.coord=[obj.var.coord;11 11 11;-11 -11 11;];                        
%                 x=obj.var.coord(:,1); y=obj.var.coord(:,2); z=obj.var.coord(:,3); plot3(x,y,z,'.');
%--------------------------------------------------------------------------
            case 15
                ver_tem = icosphere2;
                r_tem = 6.5;
                obj.var.coord = [];
                obj.var.coord(:,1) = ver_tem(:,1)*r_tem;
                obj.var.coord(:,2) = ver_tem(:,2)*r_tem;
                obj.var.coord(:,3) = ver_tem(:,3)*r_tem;
                [~,idTem]=sort(obj.var.coord(:,1));
                obj.var.coord(idTem(1:10),:)=[];
                obj.var.coord(138,:)=obj.var.coord(138,:)*1.5;
                obj.var.coord(141,:)=obj.var.coord(141,:)*1.5;
                obj.var.coord=[obj.var.coord;[-12 5 0]];
                obj.var.coord=[obj.var.coord;[-12 -5 0]];
        end
%--------------------------------------------------------------------------
           obj.var.n_coord=size(obj.var.coord,1);
        end
%==========================================================================
%==========================================================================   
        [f] = plot(obj,varargin);
end
%==========================================================================
%==========================================================================  
methods (Static)
    function mod_name = identify(varargin)
        mod_name = 'ModSubstrate';
    end
end
%==========================================================================
%==========================================================================  
end

