#include "mex.h"
#include "matrix.h"
#include "math.h"
//==========================================================================================================================
//==========================================================================================================================
double MATxVEC (int dim, double M[dim][dim], double vec[dim], int id){
    double x;
    x=0.;
    for (int i=0; i<dim; i++){
       x=x+M[id][i]*vec[i];
    }
    return x;
}
//==========================================================================================================================
//==========================================================================================================================
double potential(double r1_x, double r1_y, double r1_z, double r2_x, double r2_y, double r2_z, double k) {
   double X=r1_x-r2_x; double Y=r1_y-r2_y; double Z=r1_z-r2_z;
   return 0.5*k*(X*X+Y*Y+Z*Z);
   //return 0.5*k*(X*X*X*X+Y*Y*Y*Y+Z*Z*Z*Z);
}
double potential_LJ(double r1_x, double r1_y, double r1_z, double r2_x, double r2_y, double r2_z, double eps, int n, double r0, double dr) {
   double X=r1_x-r2_x; double Y=r1_y-r2_y; double Z=r1_z-r2_z;
   double r=sqrt(X*X+Y*Y+Z*Z);
   return eps*(pow((r0/(r-dr)),2*n)-2*pow((r0/(r-dr)),n));
   //return 0.5*k*(X*X*X*X+Y*Y*Y*Y+Z*Z*Z*Z);
}
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* coord_c_in, double* connect_in, double* pm_in,double* coord_org_in, double* a_in, double* ctrl_pt_org_in,
                double* outM, double* outM2,
                mwSize m_coord_c, mwSize m_pm, mwSize m_coord_org)
{
//----------------------------------------------------------------------------------------
double f_c[m_coord_c][6];
for (int i=0; i<m_coord_c; i++){
    f_c[i][0] = 0; f_c[i][1] = 0; f_c[i][2] = 0; f_c[i][3] = 0; f_c[i][4] = 0; f_c[i][5] = 0;
}
//----------------------------------------------------------------------------------------
double coord_c[m_coord_c][3];
for (int i=0; i<m_coord_c; i++){
    coord_c[i][0] = coord_c_in[i]; coord_c[i][1] = coord_c_in[i+m_coord_c]; coord_c[i][2] = coord_c_in[i+m_coord_c*2]; 
    //mexPrintf("%f %f %f;...\n", coord_c[i][0], coord_c[i][1], coord_c[i][2]);
}
//----------------------------------------------------------------------------------------
int connect[m_coord_c*3][2];
for (int i=0; i<m_coord_c*3; i++){
    connect[i][0] = floor(connect_in[i]+0.5)-1; connect[i][1] = floor(connect_in[i+m_coord_c*3]+0.5)-1; 
    //mexPrintf("%d %d;...\n", connect[i][0], connect[i][1]);
}
//----------------------------------------------------------------------------------------
double ctrl_pt_org[6][3];
for (int i=0; i<6; i++){
    ctrl_pt_org[i][0]=ctrl_pt_org_in[i]; ctrl_pt_org[i][1]=ctrl_pt_org_in[i+6]; ctrl_pt_org[i][2]=ctrl_pt_org_in[i+12]; 
    //mexPrintf("%f %f %f;...\n", ctrl_pt_org[i][0], ctrl_pt_org[i][1], ctrl_pt_org[i][2]);
}
//----------------------------------------------------------------------------------------
double pm[m_pm];
for (int i = 0; i < m_pm; i++) {
    pm[i] = pm_in[i];
    //mexPrintf("%f;...\n", pm[i]);
}
//----------------------------------------------------------------------------------------
double coord_org[m_coord_org][3];
for (int i=0; i<m_coord_org; i++){
    coord_org[i][0] = coord_org_in[i]; coord_org[i][1] = coord_org_in[i+m_coord_org]; coord_org[i][2] = coord_org_in[i+m_coord_org*2]; 
    //mexPrintf("%f %f %f;...\n", coord_org[i][0], coord_org[i][1], coord_org[i][2]);
}
//----------------------------------------------------------------------------------------
double a[m_coord_c][3];
for (int i=0; i<m_coord_c; i++){
    a[i][0] = a_in[i]; a[i][1] = a_in[i+m_coord_c]; a[i][2] = a_in[i+m_coord_c*2]; 
    //mexPrintf("%f %f %f;...\n", a[i][0], a[i][1], a[i][2]);
}
//========================================================================================================================== 
//========================================================================================================================== 
double O1[3][3];
double O1_tem[3][3];
double O2[3][3];
double O2_tem[3][3];
double x_tem1[3];
double a_tem1[3];
double waist1[3];
double waist1_tem[3];
double ctrlP1[3];
double ctrlP1_tem[3];
double ctrlP1_2[3];
double ctrlP1_2_tem[3];
double coord_tem1[3];
double x_tem2[3];
double a_tem2[3];
double waist2[3];
double waist2_tem[3];
double ctrlP2[3];
double ctrlP2_tem[3];
double ctrlP2_2[3];
double ctrlP2_2_tem[3];
double coord_tem2[3];
double f_tem;
int id_waist[3]; id_waist[0]=7; id_waist[1]=15; id_waist[2]=23;
double V; double V_alt;
outM2[0]=0.;
//-------------------------------------------------------------------------
int n_LJ=6; double r0_LJ=0.1; double dr_LJ=-0.1;
//-------------------------------------------------------------------------
for (int ic1=0; ic1<m_coord_c-1; ic1++){
//-------------------------------------------------------------------------
    double Phi=sqrt(a[ic1][0]*a[ic1][0]+a[ic1][1]*a[ic1][1]+a[ic1][2]*a[ic1][2]);
    double factor1=(1.-cos(Phi))/Phi/Phi;
    double factor2=1/Phi;
    O1[0][0]=a[ic1][0]*a[ic1][0]*factor1+Phi*cos(Phi)*factor2;
    O1[0][1]=a[ic1][0]*a[ic1][1]*factor1-a[ic1][2]*sin(Phi)*factor2;
    O1[0][2]=a[ic1][0]*a[ic1][2]*factor1+a[ic1][1]*sin(Phi)*factor2;
    O1[1][0]=a[ic1][1]*a[ic1][0]*factor1+a[ic1][2]*sin(Phi)*factor2;
    O1[1][1]=a[ic1][1]*a[ic1][1]*factor1+Phi*cos(Phi)*factor2;
    O1[1][2]=a[ic1][1]*a[ic1][2]*factor1-a[ic1][0]*sin(Phi)*factor2;
    O1[2][0]=a[ic1][2]*a[ic1][0]*factor1-a[ic1][1]*sin(Phi)*factor2;
    O1[2][1]=a[ic1][2]*a[ic1][1]*factor1+a[ic1][0]*sin(Phi)*factor2;
    O1[2][2]=a[ic1][2]*a[ic1][2]*factor1+Phi*cos(Phi)*factor2;
    //for (int i=0; i<3; i++){mexPrintf("%f %f %f;...\n", O[i][0], O[i][1], O[i][2]);}   
//-------------------------------------------------------------------------
    for (int ic2=ic1+1; ic2<m_coord_c; ic2++){
//-------------------------------------------------------------------------
    for (int ileg1=0; ileg1<3; ileg1++){
    if (connect[ic1*3+ileg1][0]==ic2){
    int ileg2=connect[ic1*3+ileg1][1];
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
    Phi=sqrt(a[ic2][0]*a[ic2][0]+a[ic2][1]*a[ic2][1]+a[ic2][2]*a[ic2][2]);
    factor1=(1.-cos(Phi))/Phi/Phi;
    factor2=1/Phi;
    O2[0][0]=a[ic2][0]*a[ic2][0]*factor1+Phi*cos(Phi)*factor2;
    O2[0][1]=a[ic2][0]*a[ic2][1]*factor1-a[ic2][2]*sin(Phi)*factor2;
    O2[0][2]=a[ic2][0]*a[ic2][2]*factor1+a[ic2][1]*sin(Phi)*factor2;
    O2[1][0]=a[ic2][1]*a[ic2][0]*factor1+a[ic2][2]*sin(Phi)*factor2;
    O2[1][1]=a[ic2][1]*a[ic2][1]*factor1+Phi*cos(Phi)*factor2;
    O2[1][2]=a[ic2][1]*a[ic2][2]*factor1-a[ic2][0]*sin(Phi)*factor2;
    O2[2][0]=a[ic2][2]*a[ic2][0]*factor1-a[ic2][1]*sin(Phi)*factor2;
    O2[2][1]=a[ic2][2]*a[ic2][1]*factor1+a[ic2][0]*sin(Phi)*factor2;
    O2[2][2]=a[ic2][2]*a[ic2][2]*factor1+Phi*cos(Phi)*factor2;
    //for (int i=0; i<3; i++){mexPrintf("%f %f %f;...\n", O[i][0], O[i][1], O[i][2]);}   
//-------------------------------------------------------------------------
        for (int k=0; k<3; k++) {x_tem1[k]=coord_org[id_waist[ileg1]][k];}
        for (int k=0; k<3; k++) {waist1[k]=MATxVEC(3, O1, x_tem1, k)+coord_c[ic1][k];}
        for (int k=0; k<3; k++) {x_tem1[k]=ctrl_pt_org[ileg1][k];}
        for (int k=0; k<3; k++) {ctrlP1[k]=MATxVEC(3, O1, x_tem1, k)+coord_c[ic1][k];}
        for (int k=0; k<3; k++) {x_tem1[k]=ctrl_pt_org[ileg1+3][k];}
        for (int k=0; k<3; k++) {ctrlP1_2[k]=MATxVEC(3, O1, x_tem1, k)+coord_c[ic1][k];}
        for (int k=0; k<3; k++) {x_tem2[k]=coord_org[id_waist[ileg2]][k];}
        for (int k=0; k<3; k++) {waist2[k]=MATxVEC(3, O2, x_tem2, k)+coord_c[ic2][k];}
        for (int k=0; k<3; k++) {x_tem2[k]=ctrl_pt_org[ileg2][k];}
        for (int k=0; k<3; k++) {ctrlP2[k]=MATxVEC(3, O2, x_tem2, k)+coord_c[ic2][k];}
        for (int k=0; k<3; k++) {x_tem2[k]=ctrl_pt_org[ileg2+3][k];}
        for (int k=0; k<3; k++) {ctrlP2_2[k]=MATxVEC(3, O2, x_tem2, k)+coord_c[ic2][k];}
//-------------------------------------------------------------------------       
        V= potential(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0])
          +potential(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0])
          +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
          +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
//         V= potential_LJ(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//           +potential_LJ(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//           +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
//           +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
          //mexPrintf("%f ;...\n", V);
        
        outM2[0]=outM2[0]+V;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------        
        for (int m=0; m<3; m++) {
            coord_c[ic1][m]+=pm[2];
            waist1[m]+=pm[2];
            ctrlP1[m]+=pm[2];
            ctrlP1_2[m]+=pm[2];
            V_alt= potential(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0])
                  +potential(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0])
                  +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
                  +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
//             V_alt= potential_LJ(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential_LJ(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
//                   +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
            f_tem=-(V_alt-V)/pm[2];
            f_c[ic1][m] += f_tem;            
            coord_c[ic1][m]-=pm[2];
            waist1[m]-=pm[2];
            ctrlP1[m]-=pm[2];
            ctrlP1_2[m]-=pm[2];
//-------------------------------------------------------------------------
//------------------------------------------------------------------------- 
            coord_c[ic2][m]+=pm[2];
            waist2[m]+=pm[2];
            ctrlP2[m]+=pm[2];
            ctrlP2_2[m]+=pm[2];
            V_alt= potential(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0])
                  +potential(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0])
                  +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
                  +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
//             V_alt= potential_LJ(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential_LJ(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
//                   +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
            f_tem=-(V_alt-V)/pm[2];
            f_c[ic2][m] += f_tem;            
            coord_c[ic2][m]-=pm[2];
            waist2[m]-=pm[2];
            ctrlP2[m]-=pm[2];
            ctrlP2_2[m]-=pm[2];
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------   
            a[ic1][m]+=pm[2];
            Phi=sqrt(a[ic1][0]*a[ic1][0]+a[ic1][1]*a[ic1][1]+a[ic1][2]*a[ic1][2]);
            factor1=(1.-cos(Phi))/Phi/Phi;
            factor2=1/Phi;
            O1_tem[0][0]=a[ic1][0]*a[ic1][0]*factor1+Phi*cos(Phi)*factor2;
            O1_tem[0][1]=a[ic1][0]*a[ic1][1]*factor1-a[ic1][2]*sin(Phi)*factor2;
            O1_tem[0][2]=a[ic1][0]*a[ic1][2]*factor1+a[ic1][1]*sin(Phi)*factor2;
            O1_tem[1][0]=a[ic1][1]*a[ic1][0]*factor1+a[ic1][2]*sin(Phi)*factor2;
            O1_tem[1][1]=a[ic1][1]*a[ic1][1]*factor1+Phi*cos(Phi)*factor2;
            O1_tem[1][2]=a[ic1][1]*a[ic1][2]*factor1-a[ic1][0]*sin(Phi)*factor2;
            O1_tem[2][0]=a[ic1][2]*a[ic1][0]*factor1-a[ic1][1]*sin(Phi)*factor2;
            O1_tem[2][1]=a[ic1][2]*a[ic1][1]*factor1+a[ic1][0]*sin(Phi)*factor2;
            O1_tem[2][2]=a[ic1][2]*a[ic1][2]*factor1+Phi*cos(Phi)*factor2;
            for (int k=0; k<3; k++) {x_tem1[k]=coord_org[id_waist[ileg1]][k];}
            for (int k=0; k<3; k++) {waist1_tem[k]=MATxVEC(3, O1_tem, x_tem1, k)+coord_c[ic1][k];}
            for (int k=0; k<3; k++) {x_tem1[k]=ctrl_pt_org[ileg1][k];}
            for (int k=0; k<3; k++) {ctrlP1_tem[k]=MATxVEC(3, O1_tem, x_tem1, k)+coord_c[ic1][k];}
            for (int k=0; k<3; k++) {x_tem1[k]=ctrl_pt_org[ileg1+3][k];}
            for (int k=0; k<3; k++) {ctrlP1_2_tem[k]=MATxVEC(3, O1_tem, x_tem1, k)+coord_c[ic1][k];}         
            V_alt= potential(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0])
                  +potential(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1_tem[0], waist1_tem[1], waist1_tem[2], pm[0])
                  +potential(ctrlP1_tem[0], ctrlP1_tem[1], ctrlP1_tem[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
                  +potential(ctrlP1_2_tem[0], ctrlP1_2_tem[1], ctrlP1_2_tem[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
//             V_alt= potential_LJ(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2[0], waist2[1], waist2[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential_LJ(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1_tem[0], waist1_tem[1], waist1_tem[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential(ctrlP1_tem[0], ctrlP1_tem[1], ctrlP1_tem[2], ctrlP2_2[0], ctrlP2_2[1], ctrlP2_2[2], pm[1])
//                   +potential(ctrlP1_2_tem[0], ctrlP1_2_tem[1], ctrlP1_2_tem[2], ctrlP2[0], ctrlP2[1], ctrlP2[2], pm[1]);
            f_tem=-(V_alt-V)/pm[2];
            f_c[ic1][m+3] += f_tem;   
            a[ic1][m]-=pm[2];
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
            a[ic2][m]+=pm[2];
            Phi=sqrt(a[ic2][0]*a[ic2][0]+a[ic2][1]*a[ic2][1]+a[ic2][2]*a[ic2][2]);
            factor1=(1.-cos(Phi))/Phi/Phi;
            factor2=1/Phi;
            O2_tem[0][0]=a[ic2][0]*a[ic2][0]*factor1+Phi*cos(Phi)*factor2;
            O2_tem[0][1]=a[ic2][0]*a[ic2][1]*factor1-a[ic2][2]*sin(Phi)*factor2;
            O2_tem[0][2]=a[ic2][0]*a[ic2][2]*factor1+a[ic2][1]*sin(Phi)*factor2;
            O2_tem[1][0]=a[ic2][1]*a[ic2][0]*factor1+a[ic2][2]*sin(Phi)*factor2;
            O2_tem[1][1]=a[ic2][1]*a[ic2][1]*factor1+Phi*cos(Phi)*factor2;
            O2_tem[1][2]=a[ic2][1]*a[ic2][2]*factor1-a[ic2][0]*sin(Phi)*factor2;
            O2_tem[2][0]=a[ic2][2]*a[ic2][0]*factor1-a[ic2][1]*sin(Phi)*factor2;
            O2_tem[2][1]=a[ic2][2]*a[ic2][1]*factor1+a[ic2][0]*sin(Phi)*factor2;
            O2_tem[2][2]=a[ic2][2]*a[ic2][2]*factor1+Phi*cos(Phi)*factor2;
            for (int k=0; k<3; k++) {x_tem2[k]=coord_org[id_waist[ileg2]][k];}
            for (int k=0; k<3; k++) {waist2_tem[k]=MATxVEC(3, O2_tem, x_tem2, k)+coord_c[ic2][k];}
            for (int k=0; k<3; k++) {x_tem2[k]=ctrl_pt_org[ileg2][k];}
            for (int k=0; k<3; k++) {ctrlP2_tem[k]=MATxVEC(3, O2_tem, x_tem2, k)+coord_c[ic2][k];}
            for (int k=0; k<3; k++) {x_tem2[k]=ctrl_pt_org[ileg2+3][k];}
            for (int k=0; k<3; k++) {ctrlP2_2_tem[k]=MATxVEC(3, O2_tem, x_tem2, k)+coord_c[ic2][k];}
            V_alt= potential(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2_tem[0], waist2_tem[1], waist2_tem[2], pm[0])
                  +potential(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0])
                  +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2_tem[0], ctrlP2_2_tem[1], ctrlP2_2_tem[2], pm[1])
                  +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2_tem[0], ctrlP2_tem[1], ctrlP2_tem[2], pm[1]);
//             V_alt= potential_LJ(coord_c[ic1][0], coord_c[ic1][1], coord_c[ic1][2], waist2_tem[0], waist2_tem[1], waist2_tem[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential_LJ(coord_c[ic2][0], coord_c[ic2][1], coord_c[ic2][2], waist1[0], waist1[1], waist1[2], pm[0], n_LJ, r0_LJ, dr_LJ)
//                   +potential(ctrlP1[0], ctrlP1[1], ctrlP1[2], ctrlP2_2_tem[0], ctrlP2_2_tem[1], ctrlP2_2_tem[2], pm[1])
//                   +potential(ctrlP1_2[0], ctrlP1_2[1], ctrlP1_2[2], ctrlP2_tem[0], ctrlP2_tem[1], ctrlP2_tem[2], pm[1]);
            f_tem=-(V_alt-V)/pm[2];
            f_c[ic2][m+3] += f_tem;   
            a[ic2][m]-=pm[2];
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------              
        }
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------        
    }
    }
    }
}
//========================================================================================================================== 
//========================================================================================================================== 
// for (int i=0; i<m_coord_m; i++){
//     mexPrintf("%f %f %f;...\n", f_m[i][0], f_m[i][1], f_m[i][2]);
// }
// for (int i=0; i<m_coord_c; i++){
//     mexPrintf("%f %f %f %f %f %f;...\n", f_c[i][0], f_c[i][1], f_c[i][2], f_c[i][3], f_c[i][4], f_c[i][5]);
// }

for (int i=0; i<m_coord_c; i++){
    outM[i] = f_c[i][0];          outM[i+m_coord_c] = f_c[i][1];        outM[i+m_coord_c*2] = f_c[i][2];
    outM[i+m_coord_c*3] = f_c[i][3];  outM[i+m_coord_c*4] = f_c[i][4];  outM[i+m_coord_c*5] = f_c[i][5];
}
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
    double *coord_c;         size_t m_coord_c;           coord_c = mxGetPr(prhs[0]);         m_coord_c = mxGetM(prhs[0]);
    double *connect;                                     connect = mxGetPr(prhs[1]);         
    double *pm;              size_t m_pm;                pm = mxGetPr(prhs[2]);              m_pm = mxGetM(prhs[2]);
    double *coord_org;       size_t m_coord_org;         coord_org = mxGetPr(prhs[3]);       m_coord_org=mxGetM(prhs[3]);
    double *a;                                           a = mxGetPr(prhs[4]);     
    double *ctrl_pt_org;                                 ctrl_pt_org=mxGetPr(prhs[5]);
    
                            /* create the output matrix */                                      /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_coord_c,6,mxREAL);          outMatrix = mxGetPr(plhs[0]);
    double *outMatrix2;      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);                          outMatrix2 = mxGetPr(plhs[1]);

    /* call the computational routine */
    arrayComp(coord_c,connect,pm,coord_org,a,ctrl_pt_org,
              outMatrix,outMatrix2,
              (mwSize)m_coord_c,(mwSize)m_pm,(mwSize)m_coord_org);

//     for (int i=0; i<m_J; i++){
//         mexPrintf("%f ;...\n", J[i]);
//         //mexPrintf("%f %f %f;...\n", ver[i], ver[i+m_ver], ver[i+m_ver*2]);
//     }
}
