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
//==========================================================================================================================
//==========================================================================================================================
void arrayComp (double* coord_c_in, double* coord_m_in, double* id_adp_in, double* id_rep_m_c_leg_in, double* j_T_in, double* J_in, double* pm_in, 
                double* id_td_mem_in, double* n_node_in, double* coord_org_in, double* a_in,
                double* outM, double* outM2,double* outM3,
        mwSize m_coord_c, mwSize m_coord_m, mwSize m_adp, mwSize m_rep, mwSize n_jT, mwSize m_J, mwSize m_pm, mwSize m_coord_org)
{
//----------------------------------------------------------------------------------------
double f_c[m_coord_c][6];double f_m[m_coord_m][3];
for (int i=0; i<m_coord_c; i++){
    f_c[i][0] = 0; f_c[i][1] = 0; f_c[i][2] = 0; f_c[i][3] = 0; f_c[i][4] = 0; f_c[i][5] = 0;
}
for (int i=0; i<m_coord_m; i++){
    f_m[i][0] = 0; f_m[i][1] = 0; f_m[i][2] = 0; 
}
//----------------------------------------------------------------------------------------
double coord_c[m_coord_c][3];
for (int i=0; i<m_coord_c; i++){
    coord_c[i][0] = coord_c_in[i]; coord_c[i][1] = coord_c_in[i+m_coord_c]; coord_c[i][2] = coord_c_in[i+m_coord_c*2]; 
    //mexPrintf("%f %f %f;...\n", coord_c[i][0], coord_c[i][1], coord_c[i][2]);
}
//----------------------------------------------------------------------------------------
double coord_m[m_coord_m][3];
for (int i=0; i<m_coord_m; i++){
    coord_m[i][0] = coord_m_in[i]; coord_m[i][1] = coord_m_in[i+m_coord_m]; coord_m[i][2] = coord_m_in[i+m_coord_m*2]; 
    //mexPrintf("%f %f %f;...\n", coord_m[i][0], coord_m[i][1], coord_m[i][2]);
}
//----------------------------------------------------------------------------------------
int id_adp[m_adp];
for (int i=0; i<m_adp; i++){
    id_adp[i] = floor(id_adp_in[i]+0.5)-1; 
    //mexPrintf("%d;...\n", id_adp[i]);
}
//----------------------------------------------------------------------------------------
int id_rep_m_c_leg[m_rep][3];
for (int i=0; i<m_rep; i++){
    for (int j=0; j < 3; j ++) {
    id_rep_m_c_leg[i][j] = floor(id_rep_m_c_leg_in[i+m_rep*j]+0.5)-1; 
    //mexPrintf("%d ", id_rep_m_c_leg[i][j]);
    }
    //mexPrintf(";...\n");
}
//----------------------------------------------------------------------------------------
int j_T[m_coord_m][n_jT]; int n_node[m_coord_m]; 
for (int i=0; i<m_coord_m; i++){
    n_node[i] = floor(n_node_in[i]+0.5);
    for (int j=0; j < n_jT; j ++) {
    j_T[i][j] = floor(j_T_in[i+m_coord_m*j]+0.5)-1; 
    //mexPrintf("%d ", j_T[i][j]);
    }
    //mexPrintf(";...\n");
    //mexPrintf("%d;...\n", n_node[i]);
}
//----------------------------------------------------------------------------------------
int J[m_J];
for (int i=0; i<m_J; i++){
    J[i] = floor(J_in[i]+0.5) - 1;
    //mexPrintf("%d;...\n", J[i]);
}
//----------------------------------------------------------------------------------------
double pm[m_pm];
for (int i = 0; i < m_pm; i++) {
    pm[i] = pm_in[i];
    //mexPrintf("%f;...\n", pm[i]);
}
//----------------------------------------------------------------------------------------
int id_td_mem[m_coord_c*3][m_adp];
for (int i=0; i<m_coord_c*3; i++){
    for (int j=0; j < m_adp; j++) {
    id_td_mem[i][j] = floor(id_td_mem_in[i+m_coord_c*3*j]+0.5); 
    }
}
// for (int j=0; j < m_adp; j ++) {
//     for (int i=0; i<m_coord_c*3; i++){
//         mexPrintf("%d ", id_td_mem[i][j]);
//     }
//     mexPrintf(";...\n");
// }
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
//----------------------------------------------------------------------------------------
outM3[0]=0;
//========================================================================================================================== 
//========================================================================================================================== 
double foot[3][3];
double foot_tem[3];
double O[3][3];
double x_tem[3];
double a_tem[3];
double f_tem;
double coord_tem[3];
double V; double V_alt;
double k_enhance=5.0;
for (int ic=0; ic<m_coord_c; ic++){
//-------------------------------------------------------------------------
    double Phi=sqrt(a[ic][0]*a[ic][0]+a[ic][1]*a[ic][1]+a[ic][2]*a[ic][2]);
    double factor1=(1.-cos(Phi))/Phi/Phi;
    double factor2=1/Phi;
    O[0][0]=a[ic][0]*a[ic][0]*factor1+Phi*cos(Phi)*factor2;
    O[0][1]=a[ic][0]*a[ic][1]*factor1-a[ic][2]*sin(Phi)*factor2;
    O[0][2]=a[ic][0]*a[ic][2]*factor1+a[ic][1]*sin(Phi)*factor2;
    O[1][0]=a[ic][1]*a[ic][0]*factor1+a[ic][2]*sin(Phi)*factor2;
    O[1][1]=a[ic][1]*a[ic][1]*factor1+Phi*cos(Phi)*factor2;
    O[1][2]=a[ic][1]*a[ic][2]*factor1-a[ic][0]*sin(Phi)*factor2;
    O[2][0]=a[ic][2]*a[ic][0]*factor1-a[ic][1]*sin(Phi)*factor2;
    O[2][1]=a[ic][2]*a[ic][1]*factor1+a[ic][0]*sin(Phi)*factor2;
    O[2][2]=a[ic][2]*a[ic][2]*factor1+Phi*cos(Phi)*factor2;
    //for (int i=0; i<3; i++){mexPrintf("%f %f %f;...\n", O[i][0], O[i][1], O[i][2]);}   
//-------------------------------------------------------------------------
    for (int ift=0; ift<3; ift++) {
        for (int k=0; k<3; k++) {
          x_tem[k]=coord_org[ift+48][k];
        }
        for (int k=0; k<3; k++) {
          foot[ift][k]=MATxVEC(3, O, x_tem, k)+coord_c[ic][k];
        }
    //mexPrintf("%f %f %f;...\n", foot[ift][0], foot[ift][1], foot[ift][2]);
    }
//-------------------------------------------------------------------------
    for (int j=0; j<m_adp; j++) {
        for (int ift=0; ift<3; ift++) {
            int id_tem=ic*3+ift;
            if (id_td_mem[id_tem][j]==1){
//-------------------------------------------------------------------------
                V=potential(coord_m[id_adp[j]][0], coord_m[id_adp[j]][1], coord_m[id_adp[j]][2], foot[ift][0], foot[ift][1], foot[ift][2], pm[0]);
                outM3[0]=outM3[0]+V;
//-------------------------------------------------------------------------
                for (int m=0; m<3; m++) {
                    coord_tem[m]=coord_m[id_adp[j]][m];
                }
                for (int m=0; m<3; m++) {
                    coord_tem[m]+=pm[1];
                    V_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot[ift][0], foot[ift][1], foot[ift][2], pm[0]);
                    f_tem=-(V_alt-V)/pm[1];
                    //f_tem=f_tem/((double)n_node[id_adp[j]]+1.0);
                    for (int i_node=0; i_node<n_node[id_adp[j]]; i_node++) {
                        f_m[j_T[id_adp[j]][i_node]][m] += f_tem; 
                    }
                    f_m[id_adp[j]][m] += f_tem;  
                    coord_tem[m]=coord_m[id_adp[j]][m];
                }
//-------------------------------------------------------------------------
                for (int m=0; m<3; m++) {
                    x_tem[m]=coord_org[ift+48][m];
                    foot_tem[m]=foot[ift][m];
                    a_tem[m]=a[ic][m];
                }
                for (int m=0; m<3; m++) {   
                    foot_tem[m]+=pm[1];
                    V_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot_tem[0], foot_tem[1], foot_tem[2], pm[0]);
                    f_tem=-(V_alt-V)/pm[1];
                    f_c[ic][m]+=f_tem;
                    foot_tem[m]-=pm[1];
//-------------------------------------------------------------------------
                    a_tem[m]+=pm[1];
                    //-----------------------------------------------------
                    Phi=sqrt(a_tem[0]*a_tem[0]+a_tem[1]*a_tem[1]+a_tem[2]*a_tem[2]);
                    factor1=(1.-cos(Phi))/Phi/Phi;
                    factor2=1/Phi;
                    O[0][0]=a_tem[0]*a_tem[0]*factor1+Phi*cos(Phi)*factor2;
                    O[0][1]=a_tem[0]*a_tem[1]*factor1-a_tem[2]*sin(Phi)*factor2;
                    O[0][2]=a_tem[0]*a_tem[2]*factor1+a_tem[1]*sin(Phi)*factor2;
                    O[1][0]=a_tem[1]*a_tem[0]*factor1+a_tem[2]*sin(Phi)*factor2;
                    O[1][1]=a_tem[1]*a_tem[1]*factor1+Phi*cos(Phi)*factor2;
                    O[1][2]=a_tem[1]*a_tem[2]*factor1-a_tem[0]*sin(Phi)*factor2;
                    O[2][0]=a_tem[2]*a_tem[0]*factor1-a_tem[1]*sin(Phi)*factor2;
                    O[2][1]=a_tem[2]*a_tem[1]*factor1+a_tem[0]*sin(Phi)*factor2;
                    O[2][2]=a_tem[2]*a_tem[2]*factor1+Phi*cos(Phi)*factor2;
                    //-----------------------------------------------------
                    for (int mm=0; mm<3; mm++) {
                    foot_tem[mm]=MATxVEC(3, O, x_tem, mm)+coord_c[ic][mm];
                    }
                    V_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot_tem[0], foot_tem[1], foot_tem[2], pm[0]);
                    f_tem=-(V_alt-V)/pm[1];
                    f_c[ic][m+3]+=f_tem;
                    a_tem[m]-=pm[1];
                    for (int mm=0; mm<3; mm++) {
                    foot_tem[mm]=foot[ift][mm];
                    }
                }
//-------------------------------------------------------------------------                
            }
        }
    }
//-------------------------------------------------------------------------
}
//========================================================================================================================== 
//========================================================================================================================== 
// add repulsive force
double Vrep; double Vrep_alt;
for (int i_rep=0; i_rep<m_rep; i_rep++) {
    int im=id_rep_m_c_leg[i_rep][0];
    int ic=id_rep_m_c_leg[i_rep][1];
    int ileg=id_rep_m_c_leg[i_rep][2];
    //---------------------------------------------------------------------
    double Phi=sqrt(a[ic][0]*a[ic][0]+a[ic][1]*a[ic][1]+a[ic][2]*a[ic][2]);
    double factor1=(1.-cos(Phi))/Phi/Phi;
    double factor2=1/Phi;
    O[0][0]=a[ic][0]*a[ic][0]*factor1+Phi*cos(Phi)*factor2;
    O[0][1]=a[ic][0]*a[ic][1]*factor1-a[ic][2]*sin(Phi)*factor2;
    O[0][2]=a[ic][0]*a[ic][2]*factor1+a[ic][1]*sin(Phi)*factor2;
    O[1][0]=a[ic][1]*a[ic][0]*factor1+a[ic][2]*sin(Phi)*factor2;
    O[1][1]=a[ic][1]*a[ic][1]*factor1+Phi*cos(Phi)*factor2;
    O[1][2]=a[ic][1]*a[ic][2]*factor1-a[ic][0]*sin(Phi)*factor2;
    O[2][0]=a[ic][2]*a[ic][0]*factor1-a[ic][1]*sin(Phi)*factor2;
    O[2][1]=a[ic][2]*a[ic][1]*factor1+a[ic][0]*sin(Phi)*factor2;
    O[2][2]=a[ic][2]*a[ic][2]*factor1+Phi*cos(Phi)*factor2;
    //for (int i=0; i<3; i++){mexPrintf("%f %f %f;...\n", O[i][0], O[i][1], O[i][2]);}   
//-------------------------------------------------------------------------
    for (int ift=0; ift<3; ift++) {
        for (int k=0; k<3; k++) {
          x_tem[k]=coord_org[ift+48][k];
        }
        for (int k=0; k<3; k++) {
          foot[ift][k]=MATxVEC(3, O, x_tem, k)+coord_c[ic][k];
        }
    //mexPrintf("%f %f %f;...\n", foot[ift][0], foot[ift][1], foot[ift][2]);
    }
    //---------------------------------------------------------------------
    Vrep=potential(coord_m[im][0],coord_m[im][1],coord_m[im][2],foot[ileg][0],foot[ileg][1],foot[ileg][2],pm[0]*k_enhance);
    outM3[0]=outM3[0]+Vrep;
    //---------------------------------------------------------------------
    for (int m=0; m<3; m++) {
        coord_tem[m]=coord_m[im][m];
    }
    //---------------------------------------------------------------------
    for (int m=0; m<3; m++) {
        coord_tem[m]+=pm[1];
        Vrep_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot[ileg][0], foot[ileg][1], foot[ileg][2], pm[0]*k_enhance);
        f_tem=-(Vrep_alt-Vrep)/pm[1];
        //f_tem=f_tem/((double)n_node[im]+1.0);
        for (int i_node=0; i_node<n_node[im]; i_node++) {
            f_m[j_T[im][i_node]][m] += f_tem;
        }
        f_m[im][m] += f_tem;
        coord_tem[m]-=pm[1];
    }
    //---------------------------------------------------------------------
    for (int m=0; m<3; m++) {
        x_tem[m]=coord_org[ileg+48][m];
        foot_tem[m]=foot[ileg][m];
        a_tem[m]=a[ic][m];
    }
    for (int m=0; m<3; m++) {
        foot_tem[m]+=pm[1];
        Vrep_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot_tem[0], foot_tem[1], foot_tem[2], pm[0]*k_enhance);
        f_tem=-(Vrep_alt-Vrep)/pm[1];
        f_c[ic][m]+=f_tem;
        foot_tem[m]-=pm[1];
        //-----------------------------------------------------------------
        a_tem[m]+=pm[1];
        //-----------------------------------------------------------------
        Phi=sqrt(a_tem[0]*a_tem[0]+a_tem[1]*a_tem[1]+a_tem[2]*a_tem[2]);
        factor1=(1.-cos(Phi))/Phi/Phi;
        factor2=1/Phi;
        O[0][0]=a_tem[0]*a_tem[0]*factor1+Phi*cos(Phi)*factor2;
        O[0][1]=a_tem[0]*a_tem[1]*factor1-a_tem[2]*sin(Phi)*factor2;
        O[0][2]=a_tem[0]*a_tem[2]*factor1+a_tem[1]*sin(Phi)*factor2;
        O[1][0]=a_tem[1]*a_tem[0]*factor1+a_tem[2]*sin(Phi)*factor2;
        O[1][1]=a_tem[1]*a_tem[1]*factor1+Phi*cos(Phi)*factor2;
        O[1][2]=a_tem[1]*a_tem[2]*factor1-a_tem[0]*sin(Phi)*factor2;
        O[2][0]=a_tem[2]*a_tem[0]*factor1-a_tem[1]*sin(Phi)*factor2;
        O[2][1]=a_tem[2]*a_tem[1]*factor1+a_tem[0]*sin(Phi)*factor2;
        O[2][2]=a_tem[2]*a_tem[2]*factor1+Phi*cos(Phi)*factor2;
        //-----------------------------------------------------
        for (int mm=0; mm<3; mm++) {
            foot_tem[mm]=MATxVEC(3, O, x_tem, mm)+coord_c[ic][mm];
        }
        Vrep_alt=potential(coord_tem[0], coord_tem[1], coord_tem[2], foot_tem[0], foot_tem[1], foot_tem[2], pm[0]*k_enhance);
        f_tem=-(Vrep_alt-Vrep)/pm[1];
        f_c[ic][m+3]+=f_tem;
        a_tem[m]-=pm[1];
        for (int mm=0; mm<3; mm++) {
            foot_tem[mm]=foot[ileg][mm];
        }
        //-----------------------------------------------------------------
    }
    //---------------------------------------------------------------------
}
//========================================================================================================================== 
//========================================================================================================================== 

//    for ic= 1:mod.mod{i_mod(1)}.var.n_coord
//     for j= 1:mod.mod{i_mod(2)}.var.n_adp
//         for k=1:3
//             i=(ic-1)*3+k;
//         if id_td_mem(i,j)==true
//             V(i,j)=potential(mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_adp(j),:),foot(k,:,ic),mod.TypForce.pm.(['k_' int_name]));
//         end
//         end
//     end
//    end

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
for (int i=0; i<m_coord_m; i++){
    outM2[i]=f_m[i][0]; outM2[i+m_coord_m]=f_m[i][1]; outM2[i+m_coord_m*2]=f_m[i][2]; 
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
    double *coord_m;         size_t m_coord_m;           coord_m = mxGetPr(prhs[1]);         m_coord_m = mxGetM(prhs[1]);
    double *id_adp;          size_t m_adp;               id_adp = mxGetPr(prhs[2]);          m_adp = mxGetM(prhs[2]);
    double *id_rep_m_c_leg;  size_t m_rep;               id_rep_m_c_leg = mxGetPr(prhs[3]);  m_rep = mxGetM(prhs[3]);
    double *j_T;             size_t n_jT;                j_T = mxGetPr(prhs[4]);             n_jT = mxGetN(prhs[4]);
    double *J;               size_t m_J;                 J = mxGetPr(prhs[5]);               m_J = mxGetM(prhs[5]); ///J is id_on_coord
    double *pm;              size_t m_pm;                pm = mxGetPr(prhs[6]);              m_pm = mxGetM(prhs[6]);
    double *id_td_mem;                                   id_td_mem = mxGetPr(prhs[7]);     
    double *n_node;                                      n_node = mxGetPr(prhs[8]);
    double *coord_org;       size_t m_coord_org;         coord_org = mxGetPr(prhs[9]);       m_coord_org=mxGetM(prhs[9]);
    double *a;                                           a = mxGetPr(prhs[10]);       
    
                            /* create the output matrix */                                      /* get a pointer to the real data in the output matrix */
    double *outMatrix;       plhs[0] = mxCreateDoubleMatrix((mwSize)m_coord_c,6,mxREAL);          outMatrix = mxGetPr(plhs[0]);
    double *outMatrix2;      plhs[1] = mxCreateDoubleMatrix((mwSize)m_coord_m,3,mxREAL);          outMatrix2 = mxGetPr(plhs[1]);
    double *outMatrix3;      plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);    outMatrix3 = mxGetPr(plhs[2]);

    /* call the computational routine */
    arrayComp(coord_c,coord_m,id_adp,id_rep_m_c_leg,j_T,J,pm,id_td_mem,n_node,coord_org,a,
              outMatrix,outMatrix2,outMatrix3,
              (mwSize)m_coord_c,(mwSize)m_coord_m,(mwSize)m_adp,(mwSize)m_rep,(mwSize)n_jT,(mwSize)m_J,(mwSize)m_pm,(mwSize)m_coord_org);

//     for (int i=0; i<m_J; i++){
//         mexPrintf("%f ;...\n", J[i]);
//         //mexPrintf("%f %f %f;...\n", ver[i], ver[i+m_ver], ver[i+m_ver*2]);
//     }
}
