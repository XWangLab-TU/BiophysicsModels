classdef BiophysicsApp
%==========================================================================
%--------------------------------------------------------------------------
        % BiophysicsApp collects application-based functions using
        % various modules
        %   See also model, dynamics
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------  
%==========================================================================    
    methods(Static)
        [] = TestCMEforce(mod,varargin);
        [] = tetherPull(dirRoot,dir_mod,varargin);
        [] = redBloodCell(dirRoot,dirMod,varargin);
        [dataStored] = membraneFussion(dirRoot,dirMod,varargin);
        [] = exampleFilopodia(dirRoot,dirMod,varargin);
        [] = exampleLamellipodia(dirRoot,dirMod,varargin);
        [] = exampleEndocytosis(dirRoot,dirMod,varargin);
        [] = memoryTest(dirRoot,dirMod,varargin);
        [] = ClathrinMediatedEndocytosisInit(dirRoot,dirMod,varargin);
    end
end

