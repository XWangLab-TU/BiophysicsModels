#include "mex.h"
#include "matrix.h"
#include "math.h"
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* dens, double* edg_in, double* n_node, double* pm, double* id_on_edg_in,
                double* outM, 
                mwSize m_coord, mwSize m_edg, mwSize m_pm, mwSize m_on_edg)
{
//----------------------------------------------------------------------------------------
int edg[m_edg][2];
for (int i=0; i<m_edg; i++){
    edg[i][0] = floor(edg_in[i]+0.5)-1; edg[i][1] = floor(edg_in[i+m_edg]+0.5)-1; 
}
//----------------------------------------------------------------------------------------
double k=pm[0];
int nt=floor(pm[1]+0.5)-1;
//----------------------------------------------------------------------------------------
int id_on_edg[m_on_edg];
for (int i=0; i<m_on_edg; i++){
    id_on_edg[i] = floor(id_on_edg_in[i]+0.5)-1;
}
//----------------------------------------------------------------------------------------
//========================================================================================================================== 
//========================================================================================================================== 
double dDens1[m_edg];
double dDens2[m_edg];
for (int i=0; i<nt; i++){
    for (int iedg=0; iedg<m_on_edg; iedg++){
        int j=id_on_edg[iedg];
        dDens1[j]=k*dens[edg[j][0]]/n_node[edg[j][0]];
        dDens2[j]=k*dens[edg[j][1]]/n_node[edg[j][1]];
    }
    for (int iedg=0; iedg<m_on_edg; iedg++){
        int j=id_on_edg[iedg];
        dens[edg[j][0]]=dens[edg[j][0]]-dDens1[j]+dDens2[j];
        dens[edg[j][1]]=dens[edg[j][1]]-dDens2[j]+dDens1[j];
    }
}
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
for (int i=0; i<m_coord; i++){
    outM[i] =  dens[i];
}
//----------------------------------------------------------------------------------------
}
//==========================================================================================================================
//==========================================================================================================================
/* The gateway function */
//==========================================================================================================================
//==========================================================================================================================
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
                                                /* get the value of input     get dimensions of the input matrix */
    double *dens;         size_t m_coord;           dens = mxGetPr(prhs[0]);         m_coord = mxGetM(prhs[0]); 
    double *edg;          size_t m_edg;             edg = mxGetPr(prhs[1]);          m_edg = mxGetM(prhs[1]);
    double *n_node;                                 n_node = mxGetPr(prhs[2]);          
    double *pm;           size_t m_pm;              pm = mxGetPr(prhs[3]);           m_pm=mxGetM(prhs[3]);
    double *id_on_edg;    size_t m_on_edg;          id_on_edg=mxGetPr(prhs[4]);      m_on_edg=mxGetM(prhs[4]);
    
                            /* create the output matrix */                                      /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_coord,1,mxREAL);          outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    arrayComp(dens,edg,n_node,pm,id_on_edg,
              outMatrix,
              (mwSize)m_coord,(mwSize)m_edg,(mwSize)m_pm,(mwSize)m_on_edg);

//     for (int i=0; i<m_id_on_coord; i++){
//         mexPrintf("%f ;...\n", id_on_coord[i]);
//     }
}