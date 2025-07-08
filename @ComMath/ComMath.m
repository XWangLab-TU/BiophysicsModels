classdef ComMath
    methods(Static)
        [G] = grid(x,y,z,varargin);
        [id_mesh] = getMeshID(coord,Mesh_coord, Mesh_range,Mesh_d,varargin);
        [D,id] = getMeshBasedDmin(coord1,coord2,Mesh,varargin);
        [j] = getMeshIDNeighbor(i,Mesh_range,varargin);
        [I,check]=plane_line_intersect(n,V0,P0,P1);
        Par = fit2Dcircle(XY);
        [n] = numDigitInt(N,varargin);
        xyzR = rotAxis(r,theta,phi, varargin);
        [S,c] = ellipsoid(a,c,varargin);
    end
end

