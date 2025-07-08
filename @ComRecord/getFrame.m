function [iFrames] = getFrame(dirFrames,varargin)
%--------------------------------------------------------------------------
        % getFrame returns numbering of all the '.mat' files 
        % input: 
        % dirFrames - directory of for data 
        % see variable arguments
        %   See also makeRotMov, makeStackMov, savePlot
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2024/01/25
%-------------------------------------------------------------------------- 
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dirFrames', @(x) ischar(x));
            ip.parse(dirFrames,varargin{:});
%--------------------------------------------------------------------------
        D=dir(dirFrames);
        nD=numel(D);
        iFrames=zeros(nD-2,1);
        for i=3:nD
            nChar=numel(D(i).name);
            for iChar=nChar:-1:1
                if strcmp(D(i).name(iChar),'_')==true
                    iCharStop=iChar;
                    break;
                end
            end
            iFrames(i-2)=str2double(D(i).name(iCharStop+1:end-4));
        end
        iFrames=sort(iFrames,'descend');
%--------------------------------------------------------------------------
end