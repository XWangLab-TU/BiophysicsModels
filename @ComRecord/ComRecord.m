classdef ComRecord
%==========================================================================
%--------------------------------------------------------------------------
        % ComRecord stores the data structure of directories of a serious 
        % of modules, e.g. @ModMembrane and @ ModSubstrate, for storing and
        % retrieving computational results, plotting and etc.
        %   See also model
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
%==========================================================================
    properties
       mod
       dir_all %directories
       iFrame %current number of frame for labeling images and data
       iStage %current number of stage for labeling images and data
       nRep %number of repeat for each identical-parameter simulation 
       pNameForSearch %parameter name for scan
    end
%==========================================================================    
    methods
        function obj = ComRecord(dir_data,dir_fig,dir_init,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dir_data', @(x) ischar(x));
            ip.addRequired('dir_fig', @(x) ischar(x));
            ip.addRequired('dir_init', @(x) ischar(x));
            ip.addParameter('paraAlt', struct([]), @isstruct);
            ip.addParameter('nRep', [], @isnumeric);
            ip.addParameter('mkdir_or_not', true, @islogical);
            ip.addParameter('save_rec', true, @islogical);
            ip.addParameter('deleteFiles', false, @islogical);
            ip.parse(dir_data,dir_fig,dir_init,varargin{:});
%--------------------------------------------------------------------------
            paraAlt=ip.Results.paraAlt;
            nRep=ip.Results.nRep;
            if isempty(nRep)
                spmd
                    nRep=numlabs;
                end
                nRep=numel(nRep);
            end
            obj.nRep=nRep;
            mkdir_or_not=ip.Results.mkdir_or_not;
            save_rec=ip.Results.save_rec;
            deleteFiles=ip.Results.deleteFiles;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
            obj.dir_all=cell(1,1);
            obj.dir_all{1}=struct('dir_data',dir_data,'dir_fig',dir_fig,'dir_init',dir_init);
            obj.iFrame=1;
            obj.iStage=1;
            obj.mod=cell(0,0);
            if ~isempty(paraAlt)
                obj.pNameForSearch=fieldnames(paraAlt);
            end
            obj=setParaSearch(obj,paraAlt);
            if mkdir_or_not==true
                if exist(dir_data,'dir')~=7
                    mkdir(dir_data);
                end
                if exist(dir_fig,'dir')~=7
                    mkdir(dir_fig);
                end
                if exist(dir_init,'dir')~=7
                    mkdir(dir_init);
                end
                [nX,nY]=size(obj.dir_all);
                dirName=fieldnames(obj.dir_all{1,1});
                nDirName=numel(dirName);
                idLegit=true(nDirName,1);
                for iName=1:nDirName
                    if ~strcmp(dirName{iName}(1:3),'dir')
                        idLegit(iName)=false;
                    end
                end
                dirName=dirName(idLegit);
                nDirName=numel(dirName);
                for iX=1:nX
                    for iY=1:nY
                        for iName=1:nDirName
                            if exist(obj.dir_all{iX,iY}.(dirName{iName}),'dir')~=7
                            mkdir(obj.dir_all{iX,iY}.(dirName{iName}));
                            end
                            if deleteFiles==true
                                delete([obj.dir_all{iX,iY}.(dirName{iName}) filesep '*.*']);
                            end
                        end
                    end
                end
            end
            if save_rec==true
                rec=obj;
                save([dir_data filesep 'rec_save.mat'],'rec');
            end
        end
%--------------------------------------------------------------------------        
%--------------------------------------------------------------------------
        function [] = makeMovie(obj,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj', @(x) isa(x,'ComRecord'));
            ip.addParameter('DelayTime', 0.5, @isnumeric);
            ip.addParameter('xyzLim', [], @isnumeric);
            ip.parse(obj,varargin{:});
%--------------------------------------------------------------------------
            xyzLim=ip.Results.xyzLim;
            if isempty(xyzLim)
                xyzLim=[-10 10;-10 10;-10 10];
            end
%--------------------------------------------------------------------------
            MyFolderInfo=dir(obj.dir_all.dir_fig);
            n_file=size(MyFolderInfo,1);
            id_fig=[];
            frame_seq=[];
            for i=3:n_file
                if strcmp(MyFolderInfo(i).name(end-2:end),'fig')==true
                    i_start=6;
                    i_end=size(MyFolderInfo(i).name,2)-4;
                    frame_seq=[frame_seq;str2num(MyFolderInfo(i).name(i_start:i_end))];
                    id_fig=[id_fig;i];
                end
            end
            [i_fig,id_tem]=sort(frame_seq);
            id_fig=id_fig(id_tem);
            n_fig=numel(i_fig);
            n_skip=10;
            i_skip=1;
            
            fig_name_gif=[obj.dir_all.dir_fig filesep 'mov.gif'];
            init_gif=false;
            for i = 1:n_fig
                fprintf('%d =======> %d\n', i , n_fig);  
                if i_skip==n_skip
                    i_skip=0;
                    include_frame=true;
                else
                    include_frame=false;
                end
                i_skip=i_skip+1;
                if include_frame==true
                fig_name=[obj.dir_all.dir_fig filesep MyFolderInfo(id_fig(i)).name];
                fig_name_png=[obj.dir_all.dir_fig filesep 'frame' num2str(i_fig(i)) '.png'];
                fig=openfig(fig_name);
                colorbar off
                %view([10 20]);
                xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
                set(fig,'PaperUnits','inches','PaperPosition',[0 0 20 10]);
                print(fig,fig_name_png,'-dpng','-r100');
                close(fig);
                pngImage=imread(fig_name_png);
                %convert it to an index format for the imwrite gif command
                [gifImage,gifMap] = rgb2ind(pngImage,256);
                %and finally write the gif
                %imwrite(gifImage,gifMap,['test','.gif'],'gif');
                if init_gif == false
                    imwrite(gifImage,gifMap,fig_name_gif,'gif','DelayTime',ip.Results.DelayTime, 'Loopcount',inf);
                    init_gif=true;
                else
                    imwrite(gifImage,gifMap,fig_name_gif,'gif','DelayTime',ip.Results.DelayTime,'WriteMode','append');
                end
                end
            end
        end
%==========================================================================  
        [obj] = setParaSearch(obj,paraAlt,varargin);
    end
%==========================================================================    
methods(Static)
    [] = makeRotMov(path, varargin);
    [] = makeStackMov(path, varargin);
    [] = addFrame(dirAll,M,fig,iFrame,varargin);
    [figH] = savePlot(dir_alt, fig_name,varargin);
    [iFrames] = getFrame(dirFrames,varargin);
end
end

