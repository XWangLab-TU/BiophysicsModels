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
        [V] = ModBrownianMotor_ModPolymerChain(bm,pc,varargin);
        [V,F] = Vlinear(r1,r2,pm,varargin);
        [V] = VinMembrane2hill(r1,r2,Vpm);
        [V] = VinMembrane2well(r1,r2,Vpm);
        [V,F] = VLennardJones(r1,r2,Vpm,varargin);
        [V] = VinMembrane(r,Vpm,Vcase);
        [V,F] = ModPolymerChain(r1,r2,Vpm,varargin);
        [V,F] = Vfene(r1,r2,Vpm,varargin);
        [V,F] = Vworm(cosAng,Vpm,varargin);
        [V,F] = VwormDir(dir,Vpm,varargin);
        [V] = VedgMembrane(r,Vedg);
        [V] = Vlj(r1,r2,Epw,r0,varargin);
    end
end

