#include "mex.h"
#include "matrix.h"
#include "math.h"
//==========================================================================================================================
//==========================================================================================================================
int m_w = 521288629;    /* must not be zero */
int m_z = 362436069;    /* must not be zero */

static int GetUint()
{
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return (m_z << 16) + m_w;
}

static double GetUniform()
{
    // 0 <= u < 2^32
    int u = GetUint();
    // The magic number below is 1/(2^32 + 2).
    // The result is strictly between 0 and 1.
    return (u + 1.0) * 2.328306435454494e-10+0.5;
}
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* pm_in, double* coord_in,double* id_on_coord_in,double* j_T_in,double* n_node_in,
                double* outM,
                mwSize m_pm, mwSize m_coord, mwSize m_id_on_coord, mwSize n_j_T, mwSize m_n_node)
{
//----------------------------------------------------------------------------------------
double coord[m_coord][3];
double coord_pre[m_coord][3];
for (int i=0; i<m_coord; i++){
    coord[i][0] = coord_in[i]; coord[i][1] = coord_in[i+m_coord]; coord[i][2] = coord_in[i+m_coord*2]; 
    coord_pre[i][0] = coord[i][0]; coord_pre[i][1] = coord[i][1]; coord_pre[i][2] = coord[i][2]; 
}
//----------------------------------------------------------------------------------------
int id_on_coord[m_id_on_coord];
for (int i=0; i<m_id_on_coord; i++){
    id_on_coord[i] = floor(id_on_coord_in[i]+0.5)-1; 
}
//----------------------------------------------------------------------------------------
int j_T[m_n_node][n_j_T]; int n_node[m_n_node]; 
for (int i=0; i<m_n_node; i++){
    n_node[i] = floor(n_node_in[i]+0.5);
    for (int j=0; j < n_j_T; j ++) {
    j_T[i][j] = floor(j_T_in[i+m_n_node*j]+0.5)-1; 
    }
}
//----------------------------------------------------------------------------------------
double k;
k=pm_in[0]; 
int n;
n=floor(pm_in[1]+0.5);
double l0;
l0=pm_in[2]; 
double n_on;
n_on=m_id_on_coord;
double V,Vtry;
int iVer;
double rm;
double n_perc;
double D;
//========================================================================================================================== 
//========================================================================================================================== 
for (int i=0; i<n; i++){
    rm=GetUniform();
    n_perc=rm*(n_on-1);
    int idx=floor(n_perc+0.5);
    iVer=id_on_coord[idx];
    V=0.;
    //mexPrintf("%d %d...\n", idx,iVer);
    for (int j=0; j<n_node[iVer]; j++){
        D=0.;
        for (int k=0; k<3; k++){
            D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
        }
        D=sqrt(D);
        V+=k*(D-l0)*(D-l0);
    }
    //---------------------------------------------------------------------
    bool done;
    done=false;
    for (int iAlt=0; iAlt<6; iAlt++){
        if (iAlt==0){
            coord[iVer][0]++;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if( Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][0]--;
            }
        }
        else if (iAlt==1){
            coord[iVer][0]--;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if (Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][0]++;
            }
        }
        else if (iAlt==2){
            coord[iVer][1]++;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if (Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][1]--;
            }
        }
        else if (iAlt==3){
            coord[iVer][1]--;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if (Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][1]++;
            }
        }
        else if (iAlt==4){
            coord[iVer][2]++;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if (Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][2]--;
            }
        }
        else if (iAlt==5){
            coord[iVer][2]--;
            Vtry=0.;
            for (int j=0; j<n_node[iVer]; j++){
                D=0.;
                for (int k=0; k<3; k++){
                    D+=(coord[iVer][k]-coord[j_T[iVer][j]][k])*(coord[iVer][k]-coord[j_T[iVer][j]][k]);
                }
                D=sqrt(D);
                Vtry+=k*(D-l0)*(D-l0);
            }
            if (Vtry<V) {
               V=Vtry;
               done=true;
            }
            else {
            coord[iVer][2]++;
            }
        }
        if (done==true) {break;}
    }
} //it
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
for (int i=0; i<m_coord; i++){
    outM[i] =  coord[i][0];          outM[i+m_coord] =  coord[i][1];        outM[i+m_coord*2] =  coord[i][2];
}
// mexPrintf("======%f %f %d %f;...\n", rMin, rMax, loc_relaxed,r[iEdgExo]);
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
                                                /* get the value of input       get dimensions of the input matrix */
    double *pm;            size_t m_pm;             pm = mxGetPr(prhs[0]);          m_pm=mxGetM(prhs[0]);           
    double *coord;         size_t m_coord;          coord = mxGetPr(prhs[1]);       m_coord=mxGetM(prhs[1]);
    double *id_on_coord;   size_t m_id_on_coord;    id_on_coord = mxGetPr(prhs[2]); m_id_on_coord=mxGetM(prhs[2]);
    double *j_T;           size_t n_j_T;            j_T = mxGetPr(prhs[3]);         n_j_T=mxGetN(prhs[3]);
    double *n_node;        size_t m_n_node;         n_node = mxGetPr(prhs[4]);      m_n_node=mxGetM(prhs[4]);
    
                            /* create the output matrix */                                      /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_coord,3,mxREAL);          outMatrix = mxGetPr(plhs[0]);
    //double *outMatrix2;      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);                        outMatrix2 = mxGetPr(plhs[1]);


    /* call the computational routine */
    arrayComp(pm,coord,id_on_coord,j_T,n_node,
              outMatrix,
              (mwSize)m_pm,(mwSize)m_coord,(mwSize)m_id_on_coord,(mwSize)n_j_T,(mwSize)m_n_node);

//     for (int i=0; i<m_id_on_coord; i++){
//         mexPrintf("%f ;...\n", id_on_coord[i]);
//     }
}