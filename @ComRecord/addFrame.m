function [] = addFrame(dirAll,M,fig,iFrame,varargin)
%--------------------------------------------------------------------------
        % addFrame save a new frame, as well as the corresponding matlab variable to the path set in dirAll
        % input: 
        % dirAll - directory structure for data and figure
        % M - a given @Model object
        % fig - the figure handle as the new frame
        % iFrame - the number of frame
        % optional:
        % see variable arguments
        %   See also makeRotMov, makeStackMov, savePlot
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%-------------------------------------------------------------------------- 
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dirAll', @(x) isstruct(x));
            ip.addRequired('M', @(x) isa(x,'model'));
            ip.addRequired('fig', @(x) isobject(x));
            ip.addRequired('iFrame', @(x) isnumeric(x));
            ip.addParameter('iStage', 1, @isnumeric);
            ip.addParameter('update_mod_only', true, @islogical);
            ip.addParameter('rec_fig_only', false, @islogical);
            ip.addParameter('rec_mod_only', false, @islogical);
            ip.addParameter('NmaxFrame', 1000000, @isnumeric);
            ip.addParameter('NmaxStage', 99, @isnumeric);
            ip.addParameter('PaperPosition', [0 0 5 5], @isnumeric);
            ip.parse(dirAll,M,fig,iFrame,varargin{:});
%--------------------------------------------------------------------------
            NmaxStage=ip.Results.NmaxStage;
            NmaxFrame=ip.Results.NmaxFrame;
            nDigit=ComMath.numDigitInt(int32(NmaxStage));
            paraStageFormat=['%0.' num2str(nDigit) 'd'];
            nDigit=ComMath.numDigitInt(int32(NmaxFrame));
            paraFrameFormat=['%0.' num2str(nDigit) 'd'];
            iStage=ip.Results.iStage;
            update_mod_only=ip.Results.update_mod_only;
            rec_fig_only=ip.Results.rec_fig_only;
            rec_mod_only=ip.Results.rec_mod_only;
%--------------------------------------------------------------------------            
            if iFrame>NmaxFrame*9 % 9 is fine, 10 will need one more digit
                error('iFrame > max');
            end
            if iStage>NmaxStage
                error('iStage > max');
            end
%--------------------------------------------------------------------------
            if rec_fig_only==false
            if ip.Results.rec_fig_only==false
                if ip.Results.update_mod_only==false
                    filename=[dirAll.dir_data filesep 'mod_Stg_' num2str(iStage,paraStageFormat) 'Fr_' num2str(iFrame,paraFrameFormat) '.mat'];
                else
                    filename=[dirAll.dir_data filesep 'mod_save.mat'];
                end            
                save(filename,'M','-v7.3');
            end
            end
%--------------------------------------------------------------------------
            if rec_mod_only==false
            figName=[dirAll.dir_fig filesep 'Stg_' num2str(iStage,paraStageFormat) 'Fr_' num2str(iFrame,paraFrameFormat) '.fig'];
            figNamePng=[dirAll.dir_fig filesep 'Stg_' num2str(iStage,paraStageFormat) 'Fr_' num2str(iFrame,paraFrameFormat) '.png'];
            set(fig,'PaperUnits','inches','PaperPosition',ip.Results.PaperPosition);
            print(fig,figNamePng,'-dpng','-r300');
            savefig(fig,figName);
            end
            close(fig);
%--------------------------------------------------------------------------
        end