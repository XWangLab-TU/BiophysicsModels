classdef ComPlot
    methods(Static)
        f = Sphere(X,Y,Z,i,varargin)
        f = Cylinder(x,y,z,h,r,theta,phi, varargin);
        f = PlotMod(mod, varargin);
        varargout=notBoxPlot(y,x,varargin);
        sem=notBoxPlot_SEM_calc(vect, CI);
        tint=notBoxPlot_tInterval_calc(vect, CI);
        
        PF = PassFail(fig,path,varargin);
    end
end

