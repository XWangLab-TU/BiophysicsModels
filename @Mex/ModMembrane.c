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
double Gaussian_rand () {
    double gr_1, gr_2, aux;
       do {
	gr_1 = 2. * GetUniform() - 1.;
 	gr_2 = 2. * GetUniform() - 1.;
        aux = pow(gr_1, 2.) + pow(gr_2, 2.);
      } while (aux > 1.);

      aux = sqrt(-2. * log(aux) / aux);
      gr_1 *= aux;
      gr_2 *= aux;
      return gr_1;
}
//----------------------------------------------------------------------------------------
double getMax (double* x, mwSize nx){
    double xmax;
    xmax = x[0];
    for (int i = 1; i < nx; i++) {
        if ( xmax < x[i] ) {
            xmax = x[i];
        }    
    }
    return xmax;
}
//----------------------------------------------------------------------------------------
double getMin (double* x, mwSize nx){
    double xmin;
    xmin = x[0];
    for (int i = 1; i < nx; i++) {
        if ( xmin > x[i] ) {
            xmin = x[i];
        }    
    }
    return xmin;
}
//----------------------------------------------------------------------------------------
double getMinID (double* x, mwSize nx){
    double xmin;
    int ID = 0;
    xmin = x[0];
    for (int i = 1; i < nx; i++) {
        if ( xmin > x[i] ) {
            xmin = x[i];
            ID = i;
        }    
    }
    return ID;
}
//----------------------------------------------------------------------------------------
double getCrossX (double* a, double* b){
    return a[1]*b[2]-a[2]*b[1];
}
double getCrossY (double* a, double* b){
    return -(a[0]*b[2]-a[2]*b[0]);
}
double getCrossZ (double* a, double* b){
    return a[0]*b[1]-a[1]*b[0];
}
//==========================================================================================================================
double getVtetrahedron (double* x, double* y, double* z, double* O){

double V;
double A[3]; double B[3]; double C[3];
double X[3]; double Y[3]; double Z[3];
double dirV;
double CrossRes[3];
double Edge1[3]; double Edge2[3];
 
    A[0]=x[0]; A[1]=y[0]; A[2]=z[0];
    B[0]=x[1]; B[1]=y[1]; B[2]=z[1];
    C[0]=x[2]; C[1]=y[2]; C[2]=z[2];

    X[0]=x[0]-O[0]; X[1]=x[1]-O[0]; X[2]=x[2]-O[0];
    Y[0]=y[0]-O[1]; Y[1]=y[1]-O[1]; Y[2]=y[2]-O[1];
    Z[0]=z[0]-O[2]; Z[1]=z[1]-O[2]; Z[2]=z[2]-O[2];
    
    V=fabs(-X[2]*Y[1]*Z[0]+X[1]*Y[2]*Z[0]+X[2]*Y[0]*Z[1]-X[0]*Y[2]*Z[1]-X[1]*Y[0]*Z[2]+X[0]*Y[1]*Z[2])/6;
    
    for (int iTem=0; iTem<3; iTem++) {Edge1[iTem]=B[iTem]-A[iTem]; Edge2[iTem]=B[iTem]-C[iTem];}
    
    CrossRes[0] = getCrossX(Edge1,Edge2); CrossRes[1] = getCrossY(Edge1,Edge2); CrossRes[2] = getCrossZ(Edge1,Edge2);
    
    dirV=0; for (int iTem=0; iTem<3; iTem++) { dirV=dirV+CrossRes[iTem]*(A[iTem]-O[iTem]);}
    
    if (dirV>0) {dirV=-1;}       
    else        {dirV=1;}    
    
    V=dirV*V;
    return V;
}
//==========================================================================================================================
double getForce (double r, double* r_std, int m_r_std, int* r_sub_n, int m_r_sub_n){
    int iseg=0;
    double eps = 0.0000000000000000001;
    bool found = false;
    for (int i = 0; i < m_r_sub_n-1; i++) {
        if ((r - r_std[r_sub_n[i]] > 0) && (r - r_std[r_sub_n[i+1]] < 0)) {
        iseg = i;
        found = true;
        break;
        }
    }
    if (found == false) {
    for (int i = 0; i < m_r_sub_n-1; i++) {
        if ((r - r_std[r_sub_n[i]] > -eps) && (r - r_std[r_sub_n[i+1]] < eps)) {
        iseg = i;
        found = true;
        break;
        }
    }
    }
    if (found == false) {mexPrintf("%d %f %f;\n", iseg, r, r - r_std[r_sub_n[0]]);}
    return iseg;
}
//==========================================================================================================================
double getHelfrich (int m_ver, double ver[m_ver][3], double ver_alt[m_ver][3], int* n_node, double* kH, double* A, int* idx, int m_idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
double K[m_ver][3];
double H_alt=0.; double H = 0.;
//-------------------------------------------------------------------------
double kH_alt[m_ver]; double A_alt[m_ver]; 
for (int m = 0; m < m_ver; m++) {
    kH_alt[m] = 0.; A_alt[m] = 0.;
}
//-------------------------------------------------------------------------
for (int i_ver=0; i_ver<m_idx; i_ver++){
    int i = idx[i_ver];
    H += 0.5*(2*kH[i]*2*kH[i])*A[i];
}
//-------------------------------------------------------------------------
for (int i_ver=0; i_ver<m_idx; i_ver++){
    int i = idx[i_ver];
//------------------------------------------------------------------------- 
    double c[n_node[i]][3]; double b[n_node[i]][3]; double c_d_c[n_node[i]]; double b_d_b[n_node[i]]; double b_d_c[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            c[j][k] = ver_alt[i][k]-ver_alt[j_T[i][j]][k];
            b[j][k] = ver_alt[j_T_alt[i][j]][k]-ver_alt[j_T[i][j]][k];
        } 
    }
    for (int j = 0; j < n_node[i]; j++){
        c_d_c[j] = c[j][0]*c[j][0] + c[j][1]*c[j][1] + c[j][2]*c[j][2];
        b_d_b[j] = b[j][0]*b[j][0] + b[j][1]*b[j][1] + b[j][2]*b[j][2];
        b_d_c[j] = b[j][0]*c[j][0] + b[j][1]*c[j][1] + b[j][2]*c[j][2];
    }
    //-------------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //-------------------------------------------------------------------------
    double cosine_abg[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        cosine_abg[j][0] = b_d_c[j]/(sqrt(b_d_b[j])*sqrt(c_d_c[j]));
        cosine_abg[j][1] = -b1_d_c2[j]/(sqrt(b_d_b[j])*sqrt(c2_d_c2[j]));
        cosine_abg[j][2] = c1_d_c2[j]/(sqrt(c_d_c[j])*sqrt(c2_d_c2[j]));
    }
    double id_tem_all[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            id_tem_all[j][k] = 0.5+0.5*tanh(100.*(-cosine_abg[j][k]));
        }
    }
    //-------------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
    //-------------------------------------------------------------------------
    double con_min[n_node[i]]; double id_tem[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        con_min[j] = 1;
        for (int k = 0; k < 3; k++){
            if (con_min[j] > cosine_abg[j][k]) {
                con_min[j] = cosine_abg[j][k];
            }
        }
        id_tem[j] = 0.5+0.5*tanh(100.*(con_min[j]));
    }
    //-------------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //-------------------------------------------------------------------------
    double cota[n_node[i]]; double cotb[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        cota[j] = b_d_c[j]/Dela[j];
        cotb[j] = -b1_d_c2[j]/Delb[j];
    }
    A_alt[i] = 0;
    for (int j = 0; j < n_node[i]-1; j++){
        A_alt[i] += 0.125*((cota[j]+cotb[j+1])*c_d_c[j+1])*id_tem[j]+A_T[j];
    }
    A_alt[i] += (0.125*((cota[n_node[i]-1]+cotb[0])*c_d_c[0])*id_tem[n_node[i]-1]+A_T[n_node[i]-1]);
    for (int k = 0; k < 3; k++) {
        K[i][k] = 0;
        for (int j = 0; j < n_node[i]-1; j++){
        K[i][k] += 0.5/A_alt[i]*((cota[j]+cotb[j+1])*c[j+1][k]); 
        }
        K[i][k] += 0.5/A_alt[i]*((cota[n_node[i]-1]+cotb[0])*c[0][k]); 
    }    
    kH_alt[i] = 0.5*sqrt(K[i][0]*K[i][0]+K[i][1]*K[i][1]+K[i][2]*K[i][2]);
//-------------------------------------------------------------------------
} // i
for (int i_ver=0; i_ver<m_idx; i_ver++){
    int i = idx[i_ver];
    H_alt += 0.5*(2*kH_alt[i]*2*kH_alt[i])*A_alt[i];
}
    return H_alt-H;
}
//==========================================================================================================================
double getSinglekH (int m_ver, double ver[m_ver][3], double ver_alt[m_ver][3], int* n_node, double* kH, double* A, int idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
double K[m_ver][3]; 
double A_alt[m_ver];
//-------------------------------------------------------------------------
double kH_alt; 
//-------------------------------------------------------------------------
    int i = idx;
//------------------------------------------------------------------------- 
    double c[n_node[i]][3]; double b[n_node[i]][3]; double c_d_c[n_node[i]]; double b_d_b[n_node[i]]; double b_d_c[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            c[j][k] = ver_alt[i][k]-ver_alt[j_T[i][j]][k];
            b[j][k] = ver_alt[j_T_alt[i][j]][k]-ver_alt[j_T[i][j]][k];
        } 
    }
    for (int j = 0; j < n_node[i]; j++){
        c_d_c[j] = c[j][0]*c[j][0] + c[j][1]*c[j][1] + c[j][2]*c[j][2];
        b_d_b[j] = b[j][0]*b[j][0] + b[j][1]*b[j][1] + b[j][2]*b[j][2];
        b_d_c[j] = b[j][0]*c[j][0] + b[j][1]*c[j][1] + b[j][2]*c[j][2];
    }
    //-------------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //-------------------------------------------------------------------------
    double cosine_abg[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        cosine_abg[j][0] = b_d_c[j]/(sqrt(b_d_b[j])*sqrt(c_d_c[j]));
        cosine_abg[j][1] = -b1_d_c2[j]/(sqrt(b_d_b[j])*sqrt(c2_d_c2[j]));
        cosine_abg[j][2] = c1_d_c2[j]/(sqrt(c_d_c[j])*sqrt(c2_d_c2[j]));
    }
    double id_tem_all[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            id_tem_all[j][k] = 0.5+0.5*tanh(100.*(-cosine_abg[j][k]));
        }
    }
    //-------------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
    //-------------------------------------------------------------------------
    double con_min[n_node[i]]; double id_tem[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        con_min[j] = 1;
        for (int k = 0; k < 3; k++){
            if (con_min[j] > cosine_abg[j][k]) {
                con_min[j] = cosine_abg[j][k];
            }
        }
        id_tem[j] = 0.5+0.5*tanh(100.*(con_min[j]));
    }
    //-------------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //-------------------------------------------------------------------------
    double cota[n_node[i]]; double cotb[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        cota[j] = b_d_c[j]/Dela[j];
        cotb[j] = -b1_d_c2[j]/Delb[j];
    }
    A_alt[i] = 0;
    for (int j = 0; j < n_node[i]-1; j++){
        A_alt[i] += 0.125*((cota[j]+cotb[j+1])*c_d_c[j+1])*id_tem[j]+A_T[j];
    }
    A_alt[i] += (0.125*((cota[n_node[i]-1]+cotb[0])*c_d_c[0])*id_tem[n_node[i]-1]+A_T[n_node[i]-1]);
    for (int k = 0; k < 3; k++) {
        K[i][k] = 0;
        for (int j = 0; j < n_node[i]-1; j++){
        K[i][k] += 0.5/A_alt[i]*((cota[j]+cotb[j+1])*c[j+1][k]); 
        }
        K[i][k] += 0.5/A_alt[i]*((cota[n_node[i]-1]+cotb[0])*c[0][k]); 
    }    
    kH_alt = 0.5*sqrt(K[i][0]*K[i][0]+K[i][1]*K[i][1]+K[i][2]*K[i][2]);
//-------------------------------------------------------------------------
    return kH_alt;
}
//==========================================================================================================================
double getDA (int m_ver, double ver[m_ver][3], double ver_alt[m_ver][3], int* n_node, double* A, int* idx, int m_idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
//-------------------------------------------------------------------------
double A_alt[m_ver]; 
for (int m = 0; m < m_ver; m++) {
    A_alt[m] = 0.;
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
for (int i_ver=0; i_ver<m_idx; i_ver++){
    int i = idx[i_ver];
//------------------------------------------------------------------------- 
    double c[n_node[i]][3]; double b[n_node[i]][3]; double c_d_c[n_node[i]]; double b_d_b[n_node[i]]; double b_d_c[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            c[j][k] = ver_alt[i][k]-ver_alt[j_T[i][j]][k];
            b[j][k] = ver_alt[j_T_alt[i][j]][k]-ver_alt[j_T[i][j]][k];
        } 
    }
    for (int j = 0; j < n_node[i]; j++){
        c_d_c[j] = c[j][0]*c[j][0] + c[j][1]*c[j][1] + c[j][2]*c[j][2];
        b_d_b[j] = b[j][0]*b[j][0] + b[j][1]*b[j][1] + b[j][2]*b[j][2];
        b_d_c[j] = b[j][0]*c[j][0] + b[j][1]*c[j][1] + b[j][2]*c[j][2];
    }
    //-------------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //-------------------------------------------------------------------------
    double cosine_abg[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        cosine_abg[j][0] = b_d_c[j]/(sqrt(b_d_b[j])*sqrt(c_d_c[j]));
        cosine_abg[j][1] = -b1_d_c2[j]/(sqrt(b_d_b[j])*sqrt(c2_d_c2[j]));
        cosine_abg[j][2] = c1_d_c2[j]/(sqrt(c_d_c[j])*sqrt(c2_d_c2[j]));
    }
    double id_tem_all[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            id_tem_all[j][k] = 0.5+0.5*tanh(100.*(-cosine_abg[j][k]));
        }
    }
    //-------------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
    //-------------------------------------------------------------------------
    double con_min[n_node[i]]; double id_tem[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        con_min[j] = 1;
        for (int k = 0; k < 3; k++){
            if (con_min[j] > cosine_abg[j][k]) {
                con_min[j] = cosine_abg[j][k];
            }
        }
        id_tem[j] = 0.5+0.5*tanh(100.*(con_min[j]));
    }
    //-------------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //-------------------------------------------------------------------------
    double cota[n_node[i]]; double cotb[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        cota[j] = b_d_c[j]/Dela[j];
        cotb[j] = -b1_d_c2[j]/Delb[j];
    }
    A_alt[i] = 0;
    for (int j = 0; j < n_node[i]-1; j++){
        A_alt[i] += 0.125*((cota[j]+cotb[j+1])*c_d_c[j+1])*id_tem[j]+A_T[j];
    }
    A_alt[i] += (0.125*((cota[n_node[i]-1]+cotb[0])*c_d_c[0])*id_tem[n_node[i]-1]+A_T[n_node[i]-1]);
//-------------------------------------------------------------------------
} // i

double A_new = 0.; double A_old = 0.;
for (int i_ver=0; i_ver<m_idx; i_ver++){
    int i = idx[i_ver];
    A_new += A_alt[i];
    A_old += A[i];
}
    return A_new-A_old;
}
//==========================================================================================================================
double getSingleA (int m_ver, double ver[m_ver][3], double ver_alt[m_ver][3], int* n_node, double* A, int idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
    int i = idx;
    double A_alt;
//------------------------------------------------------------------------- 
    double c[n_node[i]][3]; double b[n_node[i]][3]; double c_d_c[n_node[i]]; double b_d_b[n_node[i]]; double b_d_c[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            c[j][k] = ver_alt[i][k]-ver_alt[j_T[i][j]][k];
            b[j][k] = ver_alt[j_T_alt[i][j]][k]-ver_alt[j_T[i][j]][k];
        } 
    }
    for (int j = 0; j < n_node[i]; j++){
        c_d_c[j] = c[j][0]*c[j][0] + c[j][1]*c[j][1] + c[j][2]*c[j][2];
        b_d_b[j] = b[j][0]*b[j][0] + b[j][1]*b[j][1] + b[j][2]*b[j][2];
        b_d_c[j] = b[j][0]*c[j][0] + b[j][1]*c[j][1] + b[j][2]*c[j][2];
    }
    //-------------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //-------------------------------------------------------------------------
    double cosine_abg[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        cosine_abg[j][0] = b_d_c[j]/(sqrt(b_d_b[j])*sqrt(c_d_c[j]));
        cosine_abg[j][1] = -b1_d_c2[j]/(sqrt(b_d_b[j])*sqrt(c2_d_c2[j]));
        cosine_abg[j][2] = c1_d_c2[j]/(sqrt(c_d_c[j])*sqrt(c2_d_c2[j]));
    }
    double id_tem_all[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            id_tem_all[j][k] = 0.5+0.5*tanh(100.*(-cosine_abg[j][k]));
        }
    }
    //-------------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
    //-------------------------------------------------------------------------
    double con_min[n_node[i]]; double id_tem[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        con_min[j] = 1;
        for (int k = 0; k < 3; k++){
            if (con_min[j] > cosine_abg[j][k]) {
                con_min[j] = cosine_abg[j][k];
            }
        }
        id_tem[j] = 0.5+0.5*tanh(100.*(con_min[j]));
    }
    //-------------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //-------------------------------------------------------------------------
    double cota[n_node[i]]; double cotb[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        cota[j] = b_d_c[j]/Dela[j];
        cotb[j] = -b1_d_c2[j]/Delb[j];
    }
    A_alt = 0;
    for (int j = 0; j < n_node[i]-1; j++){
        A_alt += 0.125*((cota[j]+cotb[j+1])*c_d_c[j+1])*id_tem[j]+A_T[j];
    }
    A_alt += (0.125*((cota[n_node[i]-1]+cotb[0])*c_d_c[0])*id_tem[n_node[i]-1]+A_T[n_node[i]-1]);
    
    return A_alt;
//-------------------------------------------------------------------------
}
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* ver_in, double* pm_in, double* edg, double* face_in, double* j_T_in, double* J_in, double* n_node_in, double* T_s, double* T_e, double* dens,
        double* outM, double* outM2,double* outM3,double* outM4,double* outM5,double* outM6,
        mwSize m_ver, mwSize m_pm, mwSize m_edg, mwSize m_face, mwSize n_jT, mwSize m_J)
{
//----------------------------------------------------------------------------------------
double pm[m_pm];
//----------------------------------------------------------------------------------------
double dt_s = pm_in[1];
//----------------------------------------------------------------------------------------
double fx;double fy;double fz; double f_on_ver_x[m_ver];double f_on_ver_y[m_ver];double f_on_ver_z[m_ver];
//----------------------------------------------------------------------------------------
double ver[m_ver][3];
for (int i=0; i<m_ver; i++){
    ver[i][0] = ver_in[i]; ver[i][1] = ver_in[i+m_ver]; ver[i][2] = ver_in[i+m_ver*2]; 
    //mexPrintf("%f %f %f;...\n", ver[i][0], ver[i][1], ver[i][2]);
}
//----------------------------------------------------------------------------------------
int face[m_face][3];
for (int i=0; i<m_face; i++){
    for (int k = 0; k < 3; k ++) {
    face[i][k] = floor(face_in[i+m_face*k]+0.5) - 1;
    }
    //mexPrintf("%d %d %d;...\n", face[i][0], face[i][1], face[i][2]);
}
//----------------------------------------------------------------------------------------
int J[m_J];
for (int i=0; i<m_J; i++){
    J[i] = floor(J_in[i]+0.5) - 1;
    //mexPrintf("%d;...\n", J[i]);
}
//----------------------------------------------------------------------------------------
int j_T[m_ver][n_jT]; int j_T_alt[m_ver][n_jT]; int n_node[m_ver]; 
for (int i=0; i<m_ver; i++){
    n_node[i] = floor(n_node_in[i]+0.5);
    for (int j=0; j < n_jT; j ++) {
    j_T[i][j] = floor(j_T_in[i+m_ver*j]+0.5)-1; 
    }
}

for (int i=0; i<m_ver; i++){
    for (int j=0; j < n_jT-1; j ++) {
    j_T_alt[i][j] = j_T[i][j+1]; 
    }
    j_T_alt[i][n_node[i]-1] = j_T[i][0]; 
}
//----------------------------------------------------------------------------------------
double Vtot=0;
//========================================================================================================================== loading parameter
for (int i = 0; i < m_pm; i++) {
    pm[i] = pm_in[i];
}
//----------------------------------------------------------------------------------------
// for (int i_p = 0; i_p < m_r_sub_n; i_p++) {
//     mexPrintf("%f %f;...\n", r_std_n[r_sub_n[i_p]], f_n[r_sub_n[i_p]]);
// }
// mexPrintf("];\n");
//----------------------------------------------------------------------------------------
//==========================================================================================================================================
//================================================================================================================================== it loop
//==========================================================================================================================================
double f_b[m_ver][3]; double f_p[m_ver][3]; double f_AV[m_ver][3];
//----------------------------------------------------------------------------------------
//========================================================================================================================== loading parameter
//========================================================================================================================== force init
for (int i=0; i<m_ver; i++){
f_on_ver_x[i] = 0;
f_on_ver_y[i] = 0;
f_on_ver_z[i] = 0;
}
//========================================================================================================================== edge force
//========================================================================================================================== random force
//========================================================================================================================== external force
//========================================================================================================================== bending force
double A[m_ver]; 

double K[m_ver][3]; double kH[m_ver];
for (int m = 0; m < m_ver; m++) {
    kH[m] = 0.; A[m] = 0.;
    for (int k = 0; k < 3; k++) {
        K[m][k] = 0;
    }
}
//-------------------------------------------------------------------------
for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
//------------------------------------------------------------------------- 
    double c[n_node[i]][3]; double b[n_node[i]][3]; double c_d_c[n_node[i]]; double b_d_b[n_node[i]]; double b_d_c[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            c[j][k] = ver[i][k]-ver[j_T[i][j]][k];
            b[j][k] = ver[j_T_alt[i][j]][k]-ver[j_T[i][j]][k];
        } 
    }
    for (int j = 0; j < n_node[i]; j++){
        c_d_c[j] = c[j][0]*c[j][0] + c[j][1]*c[j][1] + c[j][2]*c[j][2];
        b_d_b[j] = b[j][0]*b[j][0] + b[j][1]*b[j][1] + b[j][2]*b[j][2];
        b_d_c[j] = b[j][0]*c[j][0] + b[j][1]*c[j][1] + b[j][2]*c[j][2];
    }
    //-------------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //-------------------------------------------------------------------------
    double cosine_abg[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        cosine_abg[j][0] = b_d_c[j]/(sqrt(b_d_b[j])*sqrt(c_d_c[j]));
        cosine_abg[j][1] = -b1_d_c2[j]/(sqrt(b_d_b[j])*sqrt(c2_d_c2[j]));
        cosine_abg[j][2] = c1_d_c2[j]/(sqrt(c_d_c[j])*sqrt(c2_d_c2[j]));
    }
    double id_tem_all[n_node[i]][3];
    for (int j = 0; j < n_node[i]; j++){
        for (int k = 0; k < 3; k++){
            id_tem_all[j][k] = 0.5+0.5*tanh(100.*(-cosine_abg[j][k]));
        }
    }
    //-------------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
//     for (int j = 0; j < n_node[i]; j++){
//     mexPrintf("%d %f %f %f;...\n", it, A_T[j], A_T_org[j], A_T_org[j]);
//     }
    //-------------------------------------------------------------------------
    double con_min[n_node[i]]; double id_tem[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        con_min[j] = 1;
        for (int k = 0; k < 3; k++){
            if (con_min[j] > cosine_abg[j][k]) {
                con_min[j] = cosine_abg[j][k];
            }
        }
        id_tem[j] = 0.5+0.5*tanh(100.*(con_min[j]));
    }
    //-------------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //-------------------------------------------------------------------------
    double cota[n_node[i]]; double cotb[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        cota[j] = b_d_c[j]/Dela[j];
        cotb[j] = -b1_d_c2[j]/Delb[j];
    }
    A[i] = 0;
    for (int j = 0; j < n_node[i]-1; j++){
        A[i] += 0.125*((cota[j]+cotb[j+1])*c_d_c[j+1])*id_tem[j]+A_T[j];
    }
    A[i] += (0.125*((cota[n_node[i]-1]+cotb[0])*c_d_c[0])*id_tem[n_node[i]-1]+A_T[n_node[i]-1]);
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < n_node[i]-1; j++){
        K[i][k] += 0.5/A[i]*((cota[j]+cotb[j+1])*c[j+1][k]); 
        }
        K[i][k] += 0.5/A[i]*((cota[n_node[i]-1]+cotb[0])*c[0][k]); 
    }    
    kH[i] = 0.5*sqrt(K[i][0]*K[i][0]+K[i][1]*K[i][1]+K[i][2]*K[i][2]);
//-------------------------------------------------------------------------
    Vtot += 0.5*pm[3]*(2*kH[i]*2*kH[i])*A[i];
//-------------------------------------------------------------------------
   // mexPrintf("%f ;...\n", A[i]);
   // mexPrintf("%f %f %f ;...\n", K[i][0],K[i][1],K[i][2]);
} // iJ
//-------------------------------------------------------------------------
double K_tem[m_ver][3]; double K_n[m_ver]; double u_K[m_ver][3];
    for (int m = 0; m < m_ver; m++) {
        K_n[m] = 0;
    }
    for (int iJ=0; iJ<m_J; iJ++){
        int m = J[iJ];
        for (int k = 0; k < 3; k++) {
            K_tem[m][k] = K[m][k]/A[m];
            K_n[m] += K_tem[m][k]*K_tem[m][k];
        }
        K_n[m] = sqrt(K_n[m]);
        for (int k = 0; k < 3; k++) {
            if (K_n[m]<0.00000000001){
                u_K[m][k] = 0.;
            }
            else{
                u_K[m][k] = K_tem[m][k]/K_n[m];
            }
        }
    }
//-------------------------------------------------------------------------
for (int iJ=0; iJ<m_J; iJ++){
    int m = J[iJ];
    int face_start1 = m; int face_start2 = j_T[m][0]; int face_start3 = j_T[m][1];
    int j_save;
    for (int j = 0; j < m_face; j ++) {
        int id_tem[3];
        for (int k = 0; k <3; k ++) {
            int diff_tem = abs(face[j][k]-face_start1);
            id_tem[k] = diff_tem;
            diff_tem = abs(face[j][k]-face_start2);
            if (diff_tem < id_tem[k]){
                id_tem[k] = diff_tem;
            }
            diff_tem = abs(face[j][k]-face_start3);
            if (diff_tem < id_tem[k]){
                id_tem[k] = diff_tem;
            }
        }
        if ((id_tem[0]+id_tem[1]+id_tem[2]) == 0) {
            j_save = j;
            break;
        }
    }
    //------------------------------
    double a1=ver[face[j_save][1]][0]-ver[face[j_save][0]][0]; double a2=ver[face[j_save][1]][1]-ver[face[j_save][0]][1]; double a3=ver[face[j_save][1]][2]-ver[face[j_save][0]][2]; 
    double b1=ver[face[j_save][2]][0]-ver[face[j_save][1]][0]; double b2=ver[face[j_save][2]][1]-ver[face[j_save][1]][1]; double b3=ver[face[j_save][2]][2]-ver[face[j_save][1]][2];
    double dir1 = a2*b3-a3*b2; double dir2 = a3*b1-a1*b3; double dir3 = a1*b2-a2*b1;
    if ((dir1*u_K[m][0]+dir2*u_K[m][1]+dir3*u_K[m][2]) < 0) {
        u_K[m][0] = -u_K[m][0]; u_K[m][1] = -u_K[m][1]; u_K[m][2] = -u_K[m][2];
    }
    //------------------------------
    //mexPrintf("%f %f %f ;...\n", u_K[m][0],u_K[m][1],u_K[m][2]);
}
    double f_px[m_ver]; double f_py[m_ver]; double f_pz[m_ver];
    for (int m = 0; m < m_ver; m++) {
        f_p[m][0] = pm[2]*u_K[m][0]*A[m]; f_p[m][1] = pm[2]*u_K[m][1]*A[m]; f_p[m][2] = pm[2]*u_K[m][2]*A[m];
    }
//-------------------------------------------------------------------------   
double ver_alt[m_ver][3];
for (int i=0; i<m_ver; i++){
    ver_alt[i][0] = ver[i][0]; ver_alt[i][1] = ver[i][1]; ver_alt[i][2] = ver[i][2]; 
}

for (int i=0; i<m_ver; i++){
    for (int k=0; k<3; k++){
        f_b[i][k] = 0.; 
    }
}
double dH;
double dr_H = 0.00001;
for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
    int m_idx = n_node[i]+1; int idx[m_idx]; 
    for (int m = 1; m < m_idx; m++) {
    idx[m] = j_T[i][m-1];
    }
    idx[0] = i;
    bool idx_keep[m_idx];
    int n_idx_keep=0;
    for (int m = 0; m < m_idx; m++) {
    idx_keep[m] = false;
        for (int jJ=0; jJ<m_J; jJ++){
          if (J[jJ]==idx[m]){
              idx_keep[m]=true;
              n_idx_keep++;
              break;
          }
        }
    }
    int idx_new[n_idx_keep]; 
    n_idx_keep=0;
    for (int m = 0; m < m_idx; m++) {
        if (idx_keep[m]==true){
            idx_new[n_idx_keep]=idx[m];
            n_idx_keep++;
        }
    }
    
    //mexPrintf("%d %d %d %d %d ;...\n", idx[0],idx[1],idx[2],idx[3],idx[4]);
    for (int k=0; k<3; k++){
     ver_alt[i][k] += dr_H; 
     dH = getHelfrich (m_ver, ver, ver_alt, n_node, kH, A, idx_new, n_idx_keep, n_jT, j_T, j_T_alt);
     f_b[i][k] += -pm[3]*dH/dr_H;
     ver_alt[i][k] -= dr_H; 
    }
    //mexPrintf("%f %f %f ;...\n", f_b[i][0],f_b[i][1],f_b[i][2]);
}

// for (int m=0; m<m_ver; m++){
//     f_on_ver_x[m] = f_on_ver_x[m] + f_b[m][0] + f_p[m][0];
//     f_on_ver_y[m] = f_on_ver_y[m] + f_b[m][1] + f_p[m][1];
//     f_on_ver_z[m] = f_on_ver_z[m] + f_b[m][2] + f_p[m][2];
//    // mexPrintf("%f %f %f  ;...\n", f_p[m][0],f_p[m][1],f_p[m][2]);
// }
//========================================================================================================================== dynamics of AV
//double k_V = 50; double k_A = 100; double k_a = 100; double V0 = 1150*0.67; double A0 = 531; double a0 = A0/floor(m_ver+0.5);
double f_AV_tem[m_ver][3];
double r_tem = 4; double k_V = pm[6]; double k_A = pm[7]; double k_a = pm[11];
double V0 = pm[8];
double A0 = pm[9]; double a0 = A0/floor(642+0.5);
for (int i=0; i<m_ver; i++){
    f_AV_tem[i][0]=0.;f_AV_tem[i][1]=0.;f_AV_tem[i][2]=0.;
    f_AV[i][0]=0.;f_AV[i][1]=0.;f_AV[i][2]=0.;
}
if ((fabs(k_V)>0.0000001) || (fabs(k_A)>0.0000001) || (fabs(k_a)>0.0000001)){
//-------------------------------------------------------------------------
double VperF[m_face];
for (int m=0; m<m_face; m++){ VperF[m] = 0.;}
double V=0;
double x[3]; double y[3]; double z[3];
double O[3]; O[0]=-100; O[1]=-100; O[2]=-100;

for (int iF=0; iF<m_face; iF++) {
    x[0]=ver[face[iF][0]][0]; x[1]=ver[face[iF][1]][0]; x[2]=ver[face[iF][2]][0];
    y[0]=ver[face[iF][0]][1]; y[1]=ver[face[iF][1]][1]; y[2]=ver[face[iF][2]][1];
    z[0]=ver[face[iF][0]][2]; z[1]=ver[face[iF][1]][2]; z[2]=ver[face[iF][2]][2];
    
    VperF[iF]=getVtetrahedron(x, y, z, O);
}
for (int iF=0; iF<m_face; iF++) {V=V+VperF[iF]; }
outM6[2]=V;
//-------------------------------------------------------------------------
double dx = 0.00001;
//-------------------------------------------------------------------------
for (int i=0; i<m_ver; i++){
    ver_alt[i][0] = ver[i][0]; ver_alt[i][1] = ver[i][1]; ver_alt[i][2] = ver[i][2];
}
double A_tot = 0.; 
double V_alt; double A_alt; double ADE_alt; 
double VperF_tem;
for (int i=0; i<m_ver; i++)
{A_tot = A_tot+A[i]; 
 Vtot += k_a*pow(A[i]/dens[i]-a0,2)/(a0);}
outM6[1]=A_tot;
//-------------------------------------------------------------------------
Vtot += k_V*pow(V-V0,2)/(V0);
Vtot += k_A*pow(A_tot-A0,2)/(A0);
//-------------------------------------------------------------------------
double alpha=0./3.1415926; double D=0.1; double dA0=4.*3.1415926*sqrt(A0/(4.*3.1415926))*D*2;
double ADE=0.;
for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
    ADE=ADE+(2*D*kH[i]*A[i]);
}
//-------------------------------------------------------------------------
r_tem=sqrt(A_tot/4/3.1415926);
// mexPrintf("%f %f %f %f %f %f;...\n", A_tot, A0, V, V0, V/(4./3.*3.1415926*pow(r_tem,3)),ADE/dA0);
for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
    for (int k = 0; k < 3; k++){
        V_alt = V;
        //-----------------------------------------------------------------
        ver_alt[i][k] += dx;
        //-----------------------------------------------------------------
        for (int iF=0; iF<m_face; iF++) {
            if ((face[iF][0]==i) || (face[iF][1]==i) || (face[iF][2]==i)) {
            x[0]=ver_alt[face[iF][0]][0]; x[1]=ver_alt[face[iF][1]][0]; x[2]=ver_alt[face[iF][2]][0];
            y[0]=ver_alt[face[iF][0]][1]; y[1]=ver_alt[face[iF][1]][1]; y[2]=ver_alt[face[iF][2]][1];
            z[0]=ver_alt[face[iF][0]][2]; z[1]=ver_alt[face[iF][1]][2]; z[2]=ver_alt[face[iF][2]][2];
            VperF_tem=getVtetrahedron(x, y, z, O);
            V_alt=V_alt-VperF[iF]+VperF_tem;
            }
        }
//         mexPrintf("%f %f;...\n", k_V,((V_alt-V0)*(V_alt-V0)/V0-(V-V0)*(V-V0)/V0)/dx);
        f_AV_tem[i][k]=-k_V*(pow(V_alt-V0,2)/(V0)-pow(V-V0,2)/(V0))/dx;
        //-----------------------------------------------------------------
        int m_idx = n_node[i]+1; int idx[m_idx]; double a_alt[m_idx]; double kH_alt[m_idx];
        for (int m = 1; m < m_idx; m++) {
            idx[m] = j_T[i][m-1];
//             mexPrintf("%d ", idx[m]);
        }
        idx[0] = i;
//         mexPrintf("%d %d;...\n", idx[0],n_node[i]);
//         double d_A = getDA (m_ver, ver, ver_alt, n_node, A, idx, m_idx, n_jT, j_T, j_T_alt); 

        double d_A=0; double dVdens=0;
        for (int m = 0; m < m_idx; m++) {
            a_alt[m]=getSingleA(m_ver, ver, ver_alt, n_node, A, idx[m], n_jT, j_T, j_T_alt);
            kH_alt[m]=getSinglekH(m_ver, ver, ver_alt, n_node, kH, A, idx[m], n_jT, j_T, j_T_alt);
            d_A=d_A+a_alt[m]-A[idx[m]];
            dVdens=dVdens+k_a*pow(a_alt[m]/dens[idx[m]]-a0,2)/(a0)-k_a*pow(A[idx[m]]/dens[idx[m]]-a0,2)/(a0);
        }
        f_AV_tem[i][k]=f_AV_tem[i][k]-dVdens/dx;
        A_alt = A_tot + d_A;
        f_AV_tem[i][k]=f_AV_tem[i][k]-k_A*(pow(A_alt-A0,2)/(A0)-pow(A_tot-A0,2)/(A0))/dx;
        //-----------------------------------------------------------------
        double d_ADE=0;
        for (int m = 0; m < m_idx; m++) {
        d_ADE=d_ADE+(2*D*kH_alt[m]*a_alt[m])-(2*D*kH[idx[m]]*A[idx[m]]);
        }

        ADE_alt=ADE+d_ADE;
        f_AV_tem[i][k]=f_AV_tem[i][k]-alpha*pm[3]*3.1415926/(2*D*D)*((ADE_alt-dA0)*(ADE_alt-dA0)/A_alt-(ADE-dA0)*(ADE-dA0)/A_tot)/dx;
        //-----------------------------------------------------------------
        ver_alt[i][k] -= dx;
        //-----------------------------------------------------------------
//         f_AV_tem[i][k]=f_AV_tem[i][k]-k_a*((A[i]+ d_A-a0)*(A[i]+ d_A-a0)/a0-(A[i]-a0)*(A[i]-a0)/a0)/dx;
        //------------------------------------------------------------------
    }
}

   //break;
//-------------------------------------------------------------------------
for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
    for (int k = 0; k < 3; k++){
//         f_AV[i][k] = f_AV_tem[i][k];
        f_AV[i][k] = f_AV[i][k]+f_AV_tem[i][k];
        double Nmean=1;
        for (int m = 0; m < n_node[i]; m++) {
            f_AV[i][k] = f_AV[i][k]+f_AV_tem[j_T[i][m]][k];
            Nmean++;
        }
        f_AV[i][k] = f_AV[i][k]/Nmean;
    }
//     mexPrintf("%f %f %f;...\n",  f_AV[i][0],f_AV[i][1],f_AV[i][2]);
}
int n_mean=(int)(pm[10]-1);
// mexPrintf("%d\n", n_mean);
for (int i_mean=0; i_mean<n_mean; i_mean++) { 
//success: i_mean<15; 
    for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
    for (int k = 0; k < 3; k++){
        f_AV_tem[i][k] = f_AV[i][k];
        f_AV[i][k]=0;
    }
    }
    for (int iJ=0; iJ<m_J; iJ++){
    int i = J[iJ];
        for (int k = 0; k < 3; k++){
            f_AV[i][k] = f_AV[i][k]+f_AV_tem[i][k];
            double Nmean=1;
            for (int m = 0; m < n_node[i]; m++) {
                f_AV[i][k] = f_AV[i][k]+f_AV_tem[j_T[i][m]][k];
                Nmean++;
            }
            f_AV[i][k] = f_AV[i][k]/Nmean;
        }
    }
}
}
//========================================================================================================================== compute dt
//==========================================================================================================================
//==========================================================================================================================  integration
//==========================================================================================================================
//================================================================================================================= it end
//==========================================================================================================================
//mexPrintf("total fext %f ;...\n",  f_ext_tot);
for (int i=0; i<m_ver; i++){
    outM[i] = f_b[i][0];  outM[i+m_ver] = f_b[i][1];  outM[i+m_ver*2] = f_b[i][2];
    outM2[i] = f_p[i][0]; outM2[i+m_ver] = f_p[i][1]; outM2[i+m_ver*2] = f_p[i][2];
    outM3[i] = u_K[i][0]; outM3[i+m_ver] = u_K[i][1]; outM3[i+m_ver*2] = u_K[i][2];
    outM4[i] = kH[i];
    outM5[i] = f_AV[i][0]; outM5[i+m_ver] = f_AV[i][1]; outM5[i+m_ver*2] = f_AV[i][2];
    //mexPrintf("%f ;...\n",  outM4[i]);
}
outM6[0]=Vtot;
//-------------------------------------------------------------------------
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
    double *ver;        size_t m_ver;           ver = mxGetPr(prhs[0]);         m_ver = mxGetM(prhs[0]);
    double *pm;         size_t m_pm;            pm = mxGetPr(prhs[1]);          m_pm = mxGetM(prhs[1]);
    double *edg;        size_t m_edg;           edg = mxGetPr(prhs[2]);         m_edg = mxGetM(prhs[2]);
    double *face;       size_t m_face;          face = mxGetPr(prhs[3]);        m_face = mxGetM(prhs[3]);
    double *j_T;        size_t n_jT;            j_T = mxGetPr(prhs[4]);         n_jT = mxGetN(prhs[4]);
    double *J;          size_t m_J;             J = mxGetPr(prhs[5]);           m_J = mxGetM(prhs[5]); ///J is id_on_coord
    double *n_node;                             n_node = mxGetPr(prhs[6]);
    double *T_s;                                T_s = mxGetPr(prhs[7]);
    double *T_e;                                T_e = mxGetPr(prhs[8]);
    double *dens;                               dens= mxGetPr(prhs[9]);
            
    
    
                            /* create the output matrix */                              /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix = mxGetPr(plhs[0]);
    double *outMatrix2;      plhs[1] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix2 = mxGetPr(plhs[1]);
    double *outMatrix3;      plhs[2] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix3 = mxGetPr(plhs[2]);
    double *outMatrix4;      plhs[3] = mxCreateDoubleMatrix((mwSize)m_ver,1,mxREAL);     outMatrix4 = mxGetPr(plhs[3]);
    double *outMatrix5;      plhs[4] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix5 = mxGetPr(plhs[4]);
    double *outMatrix6;      plhs[5] = mxCreateDoubleMatrix(3,1,mxREAL);                 outMatrix6 = mxGetPr(plhs[5]);

    /* call the computational routine */
    arrayComp(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,dens,outMatrix,outMatrix2,outMatrix3,outMatrix4,outMatrix5,outMatrix6,
              (mwSize)m_ver,(mwSize)m_pm,(mwSize)m_edg,(mwSize)m_face,(mwSize)n_jT,(mwSize)m_J);

//     for (int i=0; i<m_J; i++){
//         mexPrintf("%f ;...\n", J[i]);
//         //mexPrintf("%f %f %f;...\n", ver[i], ver[i+m_ver], ver[i+m_ver*2]);
//     }
}
