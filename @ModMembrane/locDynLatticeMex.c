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
        if ( xmax < x[i] ) {xmax = x[i];}    
    }
    return xmax;
}
//----------------------------------------------------------------------------------------
double getMin (double* x, mwSize nx){
    double xmin;
    xmin = x[0];
    for (int i = 1; i < nx; i++) {
        if ( xmin > x[i] ) {xmin = x[i];}    
    }
    return xmin;
}
//----------------------------------------------------------------------------------------
double getMinID (double* x, mwSize nx){
    double xmin;
    int ID = 0;
    xmin = x[0];
    for (int i = 1; i < nx; i++) {
        if ( xmin > x[i] ) {xmin = x[i];ID = i;}    
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
//==========================================================================================================================
double getSinglekH (int m_ver, double ver_alt[m_ver][3], int* n_node, int idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
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
double getSingleA (int m_ver, double ver_alt[m_ver][3], int* n_node, int idx, int n_jT, int j_T[m_ver][n_jT], int j_T_alt[m_ver][n_jT]) {
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
double getEedg(int m_pmVedg, double ver1x,double ver1y,double ver1z,double ver2x,double ver2y,double ver2z,
               double pmVedg[m_pmVedg]){
    double V_0=pmVedg[0];double r_1=pmVedg[1];double r_2=pmVedg[2];double rb_1=pmVedg[3];double rb_2=pmVedg[4];
    double k_w=pmVedg[5];double e_b=pmVedg[6];double e_w=pmVedg[7];double k_b=pmVedg[8];
    double r; r=sqrt((ver2x-ver1x)*(ver2x-ver1x)+(ver2y-ver1y)*(ver2y-ver1y)+(ver2z-ver1z)*(ver2z-ver1z));
    double Eedg; Eedg=V_0*(1./tanh((r-r_1)*k_w)+1./tanh((r_2-r)*k_w)-(e_b/e_w/(1+exp(-k_b*(r-rb_1)))+e_b/e_w/(1+exp(-k_b*(rb_2-r)))));
    return Eedg;
}
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* ver_in,double* pm_in,double* edg_in,double* face_in,double* j_T_in,double* J_in,double* n_node_in,double* T_s,double* T_e,double* dens,double* pmVedg,
        double* outM, double* outM2,double* outM3,double* outM4,double* outM5,double* outM6,
        mwSize m_ver, mwSize m_pm, mwSize m_edg, mwSize m_face, mwSize n_jT, mwSize m_J, mwSize m_pmVedg)
{
//----------------------------------------------------------------------------------------
double pm[m_pm];
for (int i = 0; i < m_pm; i++) {pm[i] = pm_in[i];}
//----------------------------------------------------------------------------------------
int nt=floor(pm[0]+0.5);
double dt_s=pm[1];
double P=pm[2];
double Kc=pm[3];
double Ke=pm[4];
double ddrr=pm[5];
double KV=pm[6];
double KS=pm[7];
double Vol0=pm[8];
double Surf0=pm[9];
double nAVmean=pm[10];
double Ka=pm[11];
double l0=pm[12];
double kBT=pm[13];
double relax=pm[14];
double r_1=pmVedg[1];double r_2=pmVedg[2];double rb_1=pmVedg[3];double rb_2=pmVedg[4];
// double Rremesh1=0.5*(r_1+rb_1); double Rremesh2=0.5*(r_2+rb_2); 
double Rremesh1=rb_1; double Rremesh2=rb_2; 
//----------------------------------------------------------------------------------------
double ver[m_ver][3];
for (int i=0; i<m_ver; i++){ver[i][0] = ver_in[i]; ver[i][1] = ver_in[i+m_ver]; ver[i][2] = ver_in[i+m_ver*2]; }
//----------------------------------------------------------------------------------------
int Jall[m_ver];
for (int i=0; i<m_ver; i++){Jall[i] = i;}
//----------------------------------------------------------------------------------------
int face[m_face][3];
for (int i=0; i<m_face; i++){
    for (int k = 0; k < 3; k ++) {
    face[i][k] = floor(face_in[i+m_face*k]+0.5) - 1;
    }
//     mexPrintf("%d %d %d;\n",  face[i][0],face[i][1],face[i][2]);
}
//----------------------------------------------------------------------------------------
int edg[m_edg][2];
for (int i=0; i<m_edg; i++){
    for (int k = 0; k < 2; k ++) {
    edg[i][k] = floor(edg_in[i+m_edg*k]+0.5) - 1;
    }
}
//----------------------------------------------------------------------------------------
int J[m_J];
for (int i=0; i<m_J; i++){J[i] = floor(J_in[i]+0.5) - 1;}
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
//----------------------------------------------------------------------------------------variable preparation
double Etot,Etry; double EH,EHtry,EV,EVtry,ES,EStry,EE,EEtry;
double Surf,Vol,SurfTry,VolTry,Atem,kHtem;
double VperF[m_face]; double VperF_tem;
double EEperE[m_edg];
double rm,n_perc,D;
double Rtry;
int iVer;
double A[m_ver]; double K[m_ver][3]; double kH[m_ver];
double Prob;
bool accept=false;
double n_on=m_J;
double x[3]; double y[3]; double z[3];
double O[3]; O[0]=-100; O[1]=-100; O[2]=-100;
//==========================================================================================================================================
//================================================================================================================================== it loop
//==========================================================================================================================================
bool needRemesh=false;
outM2[2]=0; //0-no remesh, 1-need remesh
for (int it=0; it<nt; it++){
for (int m = 0; m < m_ver; m++) {
    kH[m] = 0.; A[m] = 0.;
    for (int k = 0; k < 3; k++) {
        K[m][k] = 0;
    }
}
//=========================================================================
// if ((accept==true) || (it==0)){
//-------------------------------------------------------------------------initial comp of Helfrich, A
EH=0;    
for (int iJ=0; iJ<m_ver; iJ++){
    int i = Jall[iJ];
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
    //---------------------------------------------------------------------
    double b1_d_c2[n_node[i]]; double c2_d_c2[n_node[i]]; double c1_d_c2[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
        b1_d_c2[j] = b[j][0]*c[j+1][0] + b[j][1]*c[j+1][1] + b[j][2]*c[j+1][2];
        c2_d_c2[j] = c[j+1][0]*c[j+1][0] + c[j+1][1]*c[j+1][1] + c[j+1][2]*c[j+1][2];
        c1_d_c2[j] = c[j][0]*c[j+1][0] + c[j][1]*c[j+1][1] + c[j][2]*c[j+1][2];
    }
    b1_d_c2[n_node[i]-1] = b[n_node[i]-1][0]*c[0][0] + b[n_node[i]-1][1]*c[0][1] + b[n_node[i]-1][2]*c[0][2];
    c1_d_c2[n_node[i]-1] = c[n_node[i]-1][0]*c[0][0] + c[n_node[i]-1][1]*c[0][1] + c[n_node[i]-1][2]*c[0][2];
    c2_d_c2[n_node[i]-1] = c[0][0]*c[0][0] + c[0][1]*c[0][1] + c[0][2]*c[0][2];
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    double A_T_org[n_node[i]]; double A_T[n_node[i]];
    for (int j = 0; j < n_node[i]; j++){
        A_T_org[j] = 0.5*sqrt(c_d_c[j]*b_d_b[j]-b_d_c[j]*b_d_c[j]);
    }
    for (int j = 0; j < n_node[i]; j++){
        A_T[j] = 0.5*A_T_org[j]*id_tem_all[j][2];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][0];
        A_T[j] = A_T[j]+0.25*A_T_org[j]*id_tem_all[j][1];
    }
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    double Dela[n_node[i]]; double Delb[n_node[i]];
    for (int j = 0; j < n_node[i]-1; j++){
         Dela[j] = sqrt(c_d_c[j]*b_d_b[j] - b_d_c[j]*b_d_c[j]);
         Delb[j] = sqrt(c_d_c[j+1]*b_d_b[j]-b1_d_c2[j]*b1_d_c2[j]);
    }
    Dela[n_node[i]-1] = sqrt(c_d_c[n_node[i]-1]*b_d_b[n_node[i]-1] - b_d_c[n_node[i]-1]*b_d_c[n_node[i]-1]);
    Delb[n_node[i]-1] = sqrt(c_d_c[0]*b_d_b[n_node[i]-1]-b1_d_c2[n_node[i]-1]*b1_d_c2[n_node[i]-1]);
    //---------------------------------------------------------------------
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
    //---------------------------------------------------------------------
    EH += 0.5*Kc*(2*kH[i]*2*kH[i])*A[i];
}
//-------------------------------------------------------------------------initial comp of Volume
Vol=0;
for (int iF=0; iF<m_face; iF++) {
    x[0]=ver[face[iF][0]][0]; x[1]=ver[face[iF][1]][0]; x[2]=ver[face[iF][2]][0];
    y[0]=ver[face[iF][0]][1]; y[1]=ver[face[iF][1]][1]; y[2]=ver[face[iF][2]][1];
    z[0]=ver[face[iF][0]][2]; z[1]=ver[face[iF][1]][2]; z[2]=ver[face[iF][2]][2];
    VperF[iF]=getVtetrahedron(x, y, z, O);
//     mexPrintf("%d %d %d;\n",  face[iF][0],face[iF][1],face[iF][2]);
}
for (int iF=0; iF<m_face; iF++) {Vol=Vol+VperF[iF]; }
EV=KV*pow(Vol-Vol0,2)/(Vol0);
//-------------------------------------------------------------------------initial comp of Surface
Surf=0.;
for (int iJ=0; iJ<m_ver; iJ++){
    int i = Jall[iJ];
    Surf+=A[i];
}
ES=KS*pow(Surf-Surf0,2)/(Surf0);
//-------------------------------------------------------------------------initial comp of Edge
EE=0.;
for (int iEdg=0; iEdg<m_edg; iEdg++){
    EEperE[iEdg]=getEedg(m_pmVedg,ver[edg[iEdg][0]][0],ver[edg[iEdg][0]][1],ver[edg[iEdg][0]][2],ver[edg[iEdg][1]][0],ver[edg[iEdg][1]][1],ver[edg[iEdg][1]][2],
                      pmVedg);
    EE+=EEperE[iEdg];
}
//-------------------------------------------------------------------------initial comp of total    
Etot=EH+EV+ES+EE;
// }//init    
//=========================================================================loop through vertices
rm=GetUniform();
// n_perc=rm*(floor(m_ver+0.5)-1);
n_perc=rm*(n_on-1);
int idx=floor(n_perc+0.5);
// for (int idx=0;idx<m_J;idx++){
// iVer=idx;
iVer=J[idx];
int iPick=-1;
double dr=1.000000;
double force[3];
double alpha=2.;
for (int iAlt=0; iAlt<6; iAlt++){
// for (int iAlt=0; iAlt<14; iAlt++){
    //---------------------------------------------------------------------try
//     rm=GetUniform();dr=1.*rm; //dr=floor(dr+0.5);
    //--------------------------------------------------------------------- 6 nearest options
//     if (iAlt==0){ver[iVer][0]+=dr;}else if (iAlt==1){ver[iVer][1]+=dr;}else if (iAlt==2){ver[iVer][2]+=dr;}
    if (iAlt==0){ver[iVer][0]+=dr;}else if (iAlt==1){ver[iVer][0]-=dr;}else if (iAlt==2){ver[iVer][1]+=dr;}
    else if (iAlt==3){ver[iVer][1]-=dr;}else if (iAlt==4){ver[iVer][2]+=dr;}else if (iAlt==5){ver[iVer][2]-=dr;}
    //--------------------------------------------------------------------- 8 farthest options
//     else if (iAlt==6){ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;}else if (iAlt==7){ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;}
//     else if (iAlt==8){ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;}else if (iAlt==9){ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;}
//     else if (iAlt==10){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;}else if (iAlt==11){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;}
//     else if (iAlt==12){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;}else if (iAlt==13){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;}   
//     //---------------------------------------------------------------------EH,ES trial
    EHtry=EH;
    EHtry-=0.5*Kc*(2*kH[iVer]*2*kH[iVer])*A[iVer];
    SurfTry=Surf;
    SurfTry-=A[iVer];
    for (int jVer = 0; jVer < n_node[iVer]; jVer++){
        EHtry-=0.5*Kc*(2*kH[j_T[iVer][jVer]]*2*kH[j_T[iVer][jVer]])*A[j_T[iVer][jVer]];
        SurfTry-=A[j_T[iVer][jVer]];
    }
    Atem=getSingleA(m_ver, ver, n_node, iVer, n_jT, j_T, j_T_alt);
    kHtem=getSinglekH(m_ver, ver, n_node, iVer, n_jT, j_T, j_T_alt);
    EHtry+=0.5*Kc*(2*kHtem*2*kHtem)*Atem;
    SurfTry+=Atem;
    for (int jVer = 0; jVer < n_node[iVer]; jVer++){
        Atem=getSingleA(m_ver, ver, n_node, j_T[iVer][jVer], n_jT, j_T, j_T_alt);
        kHtem=getSinglekH(m_ver, ver, n_node, j_T[iVer][jVer], n_jT, j_T, j_T_alt);
        EHtry+=0.5*Kc*(2*kHtem*2*kHtem)*Atem;
        SurfTry+=Atem;
    }
    EStry= KS*pow(SurfTry-Surf0,2)/(Surf0);
    //---------------------------------------------------------------------EV trial
    EVtry=EV;
    VolTry=Vol;
    for (int iF=0; iF<m_face; iF++) {
        if ((face[iF][0]==iVer) || (face[iF][1]==iVer) || (face[iF][2]==iVer)) {
            x[0]=ver[face[iF][0]][0]; x[1]=ver[face[iF][1]][0]; x[2]=ver[face[iF][2]][0];
            y[0]=ver[face[iF][0]][1]; y[1]=ver[face[iF][1]][1]; y[2]=ver[face[iF][2]][1];
            z[0]=ver[face[iF][0]][2]; z[1]=ver[face[iF][1]][2]; z[2]=ver[face[iF][2]][2];
            VperF_tem=getVtetrahedron(x, y, z, O);
            VolTry-=VperF[iF];
            VolTry+=VperF_tem;
        }
    }
    EVtry=KV*pow(VolTry-Vol0,2)/(Vol0);
    //---------------------------------------------------------------------EE trial
    EEtry=EE;
    bool Rover=false;
    for (int iEdg=0; iEdg<m_edg; iEdg++){
        if ((edg[iEdg][0]==iVer) || (edg[iEdg][1]==iVer)) {
        EEtry-=EEperE[iEdg];
        EEtry+=getEedg(m_pmVedg,ver[edg[iEdg][0]][0],ver[edg[iEdg][0]][1],ver[edg[iEdg][0]][2],ver[edg[iEdg][1]][0],ver[edg[iEdg][1]][1],ver[edg[iEdg][1]][2],
                pmVedg);
        }
        Rtry=sqrt((ver[edg[iEdg][1]][0]-ver[edg[iEdg][0]][0])*(ver[edg[iEdg][1]][0]-ver[edg[iEdg][0]][0])+
              (ver[edg[iEdg][1]][1]-ver[edg[iEdg][0]][1])*(ver[edg[iEdg][1]][1]-ver[edg[iEdg][0]][1])+
              (ver[edg[iEdg][1]][2]-ver[edg[iEdg][0]][2])*(ver[edg[iEdg][1]][2]-ver[edg[iEdg][0]][2]));
        if ((Rtry<r_1)||(Rtry>r_2)){Rover=true;}
    }
    //---------------------------------------------------------------------Etot trial
    if (Rover==false){Etry=EHtry+EStry+EVtry+EEtry;}
    else {Etry=Etot+10000000.;}
//     force[iAlt]=-(Etry-Etot)/dr;
    accept=false;
    if (Etry<Etot) {accept=true;}
//     else {Prob=exp(-(Etry-Etot)/kBT); rm=GetUniform();if (rm<Prob){accept=true;}}
//     mexPrintf("%f %f %f %f;...\n",  EHtry,EStry,EVtry,EEtry);
//     mexPrintf("%f %f;...\n",  Etot,Etry);
//     mexPrintf("----------------------------------------------\n");
//     mexPrintf("%f %f %f %f %f ;...\n",  EH,ES,EV,EE,Etot);
//     mexPrintf("%f %f %f %f %f ;...\n",  EHtry,EStry,EVtry,EEtry,Etry);
//     mexPrintf("----------------------------------------------\n");
    //---------------------------------------------------------------------retrieve or replace
    //--------------------------------------------------------------------- nearest
//     if (iAlt==0){ver[iVer][0]-=dr;}else if (iAlt==1){ver[iVer][1]-=dr;}else if (iAlt==2){ver[iVer][2]-=dr;}
    if      (iAlt==0){ver[iVer][0]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    else if (iAlt==1){ver[iVer][0]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    else if (iAlt==2){ver[iVer][1]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    else if (iAlt==3){ver[iVer][1]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    else if (iAlt==4){ver[iVer][2]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    else if (iAlt==5){ver[iVer][2]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    //--------------------------------------------------------------------- farthest
//     else if (iAlt==6){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if (iAlt==7){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if (iAlt==8){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if (iAlt==9){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if(iAlt==10){ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if(iAlt==11){ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if(iAlt==12){ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
//     else if(iAlt==13){ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;if( accept==true) {Etot=Etry; iPick=iAlt;}}
    //---------------------------------------------------------------------
    if( accept==true){break;}
    //---------------------------------------------------------------------
} //iAlt
//     if (iPick==0)      {mexPrintf("%f %f === ",dr,ver[iVer][0]); ver[iVer][0]+=dr; mexPrintf("%f \n",ver[iVer][0]);}
    if (iPick==0)      {ver[iVer][0]+=dr;}
    else if (iPick==1) {ver[iVer][0]-=dr;}
    else if (iPick==2) {ver[iVer][1]+=dr;}
    else if (iPick==3) {ver[iVer][1]-=dr;}
    else if (iPick==4) {ver[iVer][2]+=dr;}
    else if (iPick==5) {ver[iVer][2]-=dr;}
//     else if (iPick==6) {ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;}
//     else if (iPick==7) {ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;}
//     else if (iPick==8) {ver[iVer][0]+=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;}
//     else if (iPick==9) {ver[iVer][0]+=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;}
//     else if (iPick==10){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]+=dr;}
//     else if (iPick==11){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]+=dr;}
//     else if (iPick==12){ver[iVer][0]-=dr; ver[iVer][1]+=dr; ver[iVer][2]-=dr;}
//     else if (iPick==13){ver[iVer][0]-=dr; ver[iVer][1]-=dr; ver[iVer][2]-=dr;} 
// for (int k=0;k<3;k++){
//     ver[iVer][k]+=(double)(floor(force[k]*alpha+0.5));
// }
// mexPrintf("%f %f %f %f\n",  (double)(floor(force[0]*alpha+0.5)),(double)(floor(force[1]*alpha+0.5)),(double)(floor(force[2]*alpha+0.5)),Etot);
// }//idx
//-------------------------------------------------------------------------if need remesh 
needRemesh=false;
// mexPrintf("%f;...\n",  relax);
// if (relax<0) {
// double RtryMin=10000.,RtryMax=0.;
// for (int iEdg=0; iEdg<m_edg; iEdg++){
//     Rtry=sqrt((ver[edg[iEdg][1]][0]-ver[edg[iEdg][0]][0])*(ver[edg[iEdg][1]][0]-ver[edg[iEdg][0]][0])+
//               (ver[edg[iEdg][1]][1]-ver[edg[iEdg][0]][1])*(ver[edg[iEdg][1]][1]-ver[edg[iEdg][0]][1])+
//               (ver[edg[iEdg][1]][2]-ver[edg[iEdg][0]][2])*(ver[edg[iEdg][1]][2]-ver[edg[iEdg][0]][2]));
//     if (Rtry>RtryMax){RtryMax=Rtry;}
//     if (Rtry<RtryMin){RtryMin=Rtry;}
//     if (((Rtry<Rremesh1) || (Rtry>Rremesh2)) && (it>-1)) {needRemesh=true; break;};
// }
// }
outM2[1]=it+1;
// if (needRemesh==true) {mexPrintf("need remesh\n"); outM2[2]=1; break;}
//=========================================================================
} //it
mexPrintf("remesh done====================================================\n");
//==========================================================================================================================
//================================================================================================================= it end
//==========================================================================================================================
//=========================================================================output
//mexPrintf("total fext %f ;...\n",  f_ext_tot);
for (int i=0; i<m_ver; i++){
    outM[i] = ver[i][0];  outM[i+m_ver] = ver[i][1];  outM[i+m_ver*2] = ver[i][2];
//     outM2[i] = f_p[i][0]; outM2[i+m_ver] = f_p[i][1]; outM2[i+m_ver*2] = f_p[i][2];
//     outM3[i] = u_K[i][0]; outM3[i+m_ver] = u_K[i][1]; outM3[i+m_ver*2] = u_K[i][2];
//     outM4[i] = kH[i];
//     outM5[i] = f_AV[i][0]; outM5[i+m_ver] = f_AV[i][1]; outM5[i+m_ver*2] = f_AV[i][2];
    //mexPrintf("%f ;...\n",  outM4[i]);
}
outM2[0]=Etot;
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
//---------------------------------------------------------------------------------------------------------------------------------------------------    
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
    double *pmVedg;     size_t m_pmVedg;        pmVedg=mxGetPr(prhs[10]);       m_pmVedg=mxGetM(prhs[10]);
//---------------------------------------------------------------------------------------------------------------------------------------------------    
                            /* create the output matrix */                              /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix = mxGetPr(plhs[0]);
    double *outMatrix2;      plhs[1] = mxCreateDoubleMatrix(3,1,mxREAL);                 outMatrix2 = mxGetPr(plhs[1]); /* Etot, it, remesh?*/ 
    double *outMatrix3;      plhs[2] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix3 = mxGetPr(plhs[2]);
    double *outMatrix4;      plhs[3] = mxCreateDoubleMatrix((mwSize)m_ver,1,mxREAL);     outMatrix4 = mxGetPr(plhs[3]);
    double *outMatrix5;      plhs[4] = mxCreateDoubleMatrix((mwSize)m_ver,3,mxREAL);     outMatrix5 = mxGetPr(plhs[4]);
    double *outMatrix6;      plhs[5] = mxCreateDoubleMatrix(3,1,mxREAL);                 outMatrix6 = mxGetPr(plhs[5]);
//---------------------------------------------------------------------------------------------------------------------------------------------------
    /* call the computational routine */
    arrayComp(ver,pm,edg,face,j_T,J,n_node,T_s,T_e,dens,pmVedg,outMatrix,outMatrix2,outMatrix3,outMatrix4,outMatrix5,outMatrix6,
              (mwSize)m_ver,(mwSize)m_pm,(mwSize)m_edg,(mwSize)m_face,(mwSize)n_jT,(mwSize)m_J,(mwSize)m_pmVedg);
}