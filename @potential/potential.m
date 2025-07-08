classdef potential
%==========================================================================
%--------------------------------------------------------------------------
        % potential collects various potential functions
        %   See also model, dynamics
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/03
%--------------------------------------------------------------------------  
%==========================================================================    
    methods(Static)
        [V,F] = Vlinear(r1,r2,pm,varargin);
        [V] = Vspring(r1,r2,l0,k,varargin);
        [V] = Vising(spin1,spin2,J,varargin);
        [V] = VLennard(r1,r2,Epw,r0,varargin);
    end
end

