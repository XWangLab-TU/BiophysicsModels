classdef ComExternal
    methods(Static)
        %% remesh
        [vnew, fnew]=remesh_cleanpatch(V, F);
        [vnew, fnew,temp] = remesh_edgecollaps( vnew, fnew, sizzz,voriginal,foriginal );
        [projections]=remesh_project(vS,fS,vT,fT);
        [vnew, fnew, meanedge, stdev]=remesh_remesher(V, F, edgelength, iterations);
        [vnew, fnew, temp] = remesh_removebadtriangles( vnew, fnew ,voriginal,foriginal);
        [vnew,fnew] = remesh_subdividelarge( vnew, fnew,flag,voriginal,foriginal );
        [mfRefinedMesh, mnTriangulation] = remesh_LoopSubdivisionLimited(mfMeshPoints, mnTriangulation, fMinResolution, vbBoundaryEdges);
    end
end