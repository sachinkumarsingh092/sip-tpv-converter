#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <wcslib/dis.h>
#include <wcslib/lin.h>
#include <wcslib/wcslib.h>

#include <gsl/gsl_linalg.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/permutation.h>





/***********************************/

#define GAL_WCS_MAX_PVSIZE    40
#define GAL_WCS_MAX_POLYORDER 8

#define max(a,b)                \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b)                \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })




enum distortion_type {
  DIS_TPV = 1,
  DIS_SIP = 2,
  DIS_DSS = 3,
  DIS_WAT = 4,
  DIS_TPD = 5
};


/***********************************/







/* add filename etc. parameters later.*/
static void
gal_wcs_get_params(struct wcsprm *wcs, double cd[2][2], \
                   double *pv1, double *pv2)
{
  size_t i,j,k,index=0;
  size_t pv_m=0;
  // double cd[2][2]={0};
  // double pv1[GAL_WCS_MAX_PVSIZE] = {0};
  // double pv2[GAL_WCS_MAX_PVSIZE] = {0};
  // printf("==>%lf\n",(wcs->cd)[0]);

  /* TODO:
     [ ]Throw error if a valuse of fits file is not present. 
     [x]Also account for the missing PV terms in final arrays.
     */

  for(i=0,k=0; i<2; i++)
    for(j=0; j<2; j++)
      {
        /* If a value is present store it, else store 0.*/
        if((wcs->cd)[i] != 0)   cd[i][j]=(wcs->cd)[k++];
        else                    cd[i][j]=0;

        /* For a check:
        printf("CD%ld_%ld\t%.10lf\n", i+1, j+1, cd[i][j]);
        */
      }

  for(j=0; j<GAL_WCS_MAX_PVSIZE; j++)
  {
    if(wcs->pv[j].i == 1)
      {
        /*pv_m is used to check the index in the header.*/
        pv_m=wcs->pv[j].m;

        /* `index` is the index of the pv* array.*/
        index = pv_m;

        // printf("%d %d\n", index, pv_m);
        if( wcs->pv[pv_m].value != 0 && index == pv_m ) 
          pv1[index]=wcs->pv[j].value;
        else
          pv1[index]=0;

        /* For a check
        printf("PV1_%d\t%.10f\n", index, pv1[index]);
        */
      }
    else if(wcs->pv[j].i == 2)
      {
        /*pv_m is used to check the index in the header.*/
        pv_m=wcs->pv[j].m;

        /* `index` is the index of the pv* array.*/
        index = pv_m;
        // printf("%d %d\n", index, pv_m);
        if( wcs->pv[pv_m].value != 0 && index == pv_m ) 
          pv2[index]=wcs->pv[j].value;
        else                    
          pv2[index]=0;

        /* For a check
        printf("PV2_%d\t%.10f\n", pv_m, pv2[k]);
        */
      }
    else
      break;
  }

  /* To check a full array:
  for(i = 0; i < 40; i++)
    printf("PV2_%ld\t%.10f\n", i, pv2[i]);
  */

}






static void
gal_wcs_uv_sipeq(struct wcsprm *wcs)
{
  /*TODO*/
}






static void
gal_wcs_compute_tpv(double *pv1, double *pv2, \
                    double tpvx[8][8], double tpvy[8][8])
{
  // double tpvx[8][8]={0};
  // double tpvy[8][8]={0};
  size_t i, j;


  /* The indexes of the tpv* matrices are the exponents of x and y in 
     a TPV distortion equation leaving out radial terms 
     PV[1,3], PV[1,11], PV[1,23], PV[1,39] */

  /****** tpvx matrix *****/
  /* Order 2 terms.*/
  tpvx[0][0] = pv1[0];
  tpvx[1][0] = pv1[1];
  tpvx[0][1] = pv1[2];
  tpvx[2][0] = pv1[4];
  tpvx[1][1] = pv1[5];
  tpvx[0][2] = pv1[6];

  /*Order 3 terms.*/
  tpvx[3][0] = pv1[7];
  tpvx[2][1] = pv1[8];
  tpvx[1][2] = pv1[9];
  tpvx[0][3] = pv1[10];

  /*Order 4 terms.*/
  tpvx[4][0] = pv1[12];
  tpvx[3][1] = pv1[13];
  tpvx[2][2] = pv1[14];
  tpvx[1][3] = pv1[15];
  tpvx[0][4] = pv1[16];

  /*Order 5 terms.*/
  tpvx[5][0] = pv1[17];
  tpvx[4][1] = pv1[18];
  tpvx[3][2] = pv1[19];
  tpvx[2][3] = pv1[20];
  tpvx[1][4] = pv1[21];
  tpvx[0][5] = pv1[22];

  /*Order 6 terms.*/
  tpvx[6][0] = pv1[24];
  tpvx[5][1] = pv1[25];
  tpvx[4][2] = pv1[26];
  tpvx[3][3] = pv1[27];
  tpvx[2][4] = pv1[28];
  tpvx[1][5] = pv1[29];
  tpvx[0][6] = pv1[30];

  /*Order 7 terms.*/
  tpvx[7][0] = pv1[31];
  tpvx[6][1] = pv1[32];
  tpvx[5][2] = pv1[33];
  tpvx[4][3] = pv1[34];
  tpvx[3][4] = pv1[35];
  tpvx[2][5] = pv1[36];
  tpvx[1][6] = pv1[37];
  tpvx[0][7] = pv1[38];



  /* fromed by x <-> y conversion of above matrix*/
  /***** tpvy matrix ******/
  /* Order 2 terms.*/
  tpvy[0][1] = pv2[1];
  tpvy[0][0] = pv2[0];
  tpvy[1][0] = pv2[2];
  tpvy[0][2] = pv2[4];
  tpvy[1][1] = pv2[5];
  tpvy[2][0] = pv2[6];

  /*Order 3 terms.*/
  tpvy[0][3] = pv2[7];
  tpvy[1][2] = pv2[8];
  tpvy[2][1] = pv2[9];
  tpvy[3][0] = pv2[10];

  /*Order 4 terms.*/
  tpvy[0][4] = pv2[12];
  tpvy[1][3] = pv2[13];
  tpvy[2][2] = pv2[14];
  tpvy[3][1] = pv2[15];
  tpvy[4][0] = pv2[16];

  /*Order 5 terms.*/
  tpvy[0][5] = pv2[17];
  tpvy[1][4] = pv2[18];
  tpvy[2][3] = pv2[19];
  tpvy[3][2] = pv2[20];
  tpvy[4][1] = pv2[21];
  tpvy[5][0] = pv2[22];

  /*Order 6 terms.*/
  tpvy[0][6] = pv2[24];
  tpvy[1][5] = pv2[25];
  tpvy[2][4] = pv2[26];
  tpvy[3][3] = pv2[27];
  tpvy[4][2] = pv2[28];
  tpvy[5][1] = pv2[29];
  tpvy[6][0] = pv2[30];

  /*Order 7 terms.*/
  tpvy[0][7] = pv2[31];
  tpvy[1][6] = pv2[32];
  tpvy[2][5] = pv2[33];
  tpvy[3][4] = pv2[34];
  tpvy[4][3] = pv2[35];
  tpvy[5][2] = pv2[36];
  tpvy[6][1] = pv2[37];
  tpvy[7][0] = pv2[38];

  /* For a check:
  for(i=0; i<8; i++)
    for(j=0;j<8;j++)
      {
        printf("tpvx%ld_%ld \t %.10f\n",   i, j, tpvx[i][j]);
        printf("tpvy%ld_%ld \t %.10f\n\n", i, j, tpvy[i][j]);
      }
  */

}






static void
gal_wcs_real_tpveq(double cd[2][2], double tpvu[8][8], double tpvv[8][8])
{
  /* tpvu and tpvv are u-v distortion coefficientsin TPV convention. */
  int a;
  size_t i,j;
  double determinant=0;
  struct wcsprm *wcs=NULL;
  double cd_inv[2][2]={0};
  double tpv1[8][8]={0}, tpv2[8][8]={0};
  double pv1[GAL_WCS_MAX_PVSIZE] = {0};
  double pv2[GAL_WCS_MAX_PVSIZE] = {0};

  wcs=gal_wcs_read("test-pv.fits", "1", 0, 0, &a); 

  gal_wcs_get_params(wcs, cd, pv1, pv2);

  gal_wcs_compute_tpv(pv1, pv2, tpv1, tpv2);

  /* We will find matrix tpvu and tpvv by finding inverse of
     CD matrix and multiplying with tpv* matrix. 
     For inverse of a 2x2 matrix we use the below trasformations:
              inverse(|a  b|) =  1 *|d  -b| 
                      |c  d|    |A| |-c  a|
      where |A| is the determinant of the matrix which is calculated by:
                          |A| = a*d-b*c.
    */
  determinant = cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0];
  /* For a check:
  printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", cd[0][0], cd[0][1], \
                                             cd[1][0], cd[1][1]);
  printf("%.10lf\n", determinant);
  */

  /* Inverse matrix */
  cd_inv[0][0]=cd[1][1]/determinant;    /* a */
  cd_inv[0][1]=-1*cd[0][1]/determinant; /* b */
  cd_inv[1][0]=-1*cd[1][0]/determinant; /* c */
  cd_inv[1][1]=cd[0][0]/determinant;    /* d */

  /* For matrix tpvv and tpvu, we have to use the following 
     matrix equation:

                  |tpvu| = cd_inv*|tpv1|
                  |tpvv|          |tpv2|
      */

  for(i=0; i<8; i++){
    for(j=0; j<8; j++){
      tpvu[i][j]=cd_inv[0][0]*tpv1[i][j]+cd_inv[0][1]*tpv2[i][j];
      tpvv[i][j]=cd_inv[1][0]*tpv1[i][j]+cd_inv[1][1]*tpv2[i][j];

      /*For a check:
      printf("%.10lf, %.10lf\n", tpvu[i][j], tpvv[i][j]);
      */
    }
  }

}





/* Calculate the SIP coefficients from CD matrix
    parameters and PV coefficients. */
static double
gal_wcs_calcsip(size_t axis, size_t m, size_t n, \
                double tpvu[8][8], double tpvv[8][8])
{
  double sip_coeff;
  if(axis == 1)        sip_coeff=tpvu[m][n];
  else if(axis == 2)   sip_coeff=tpvv[m][n];
  else
    error(EXIT_FAILURE, 0, "%s: axis does not exists! ",
          __func__);

  if( (axis == 1) && (m == 1) && (n == 0) )
        sip_coeff = sip_coeff - 1.0;
  else if( (axis == 2) && (m == 0) && (n == 1) )
        sip_coeff = sip_coeff - 1.0;
    
  return sip_coeff;

}




void
gal_wcs_add_sipkeywords(double tpvu[8][8], double tpvv[8][8])
{
  size_t i,j;
  double val=0;
  size_t a_order=0, b_order=0;

  for(i=0; i<8; i++)
    for(j=0; j<8; j++)
      {
        /*For axis = 1*/
        val=gal_wcs_calcsip(1, i, j, tpvu, tpvv);
        if(val != 0)
          {
            a_order=max(a_order, max(i,j));
          }
        
        /*For axis = 2*/
        val=gal_wcs_calcsip(2, i, j, tpvu, tpvv);
        if(val != 0)
          {
            b_order=max(b_order, max(i,j));
          }
      }

}






// int
// gal_wcs_distortion_read(struct wcsprm *wcs, double *distortion)
// {
// }






int main(){
    struct wcsprm *wcs;
    int a;
    // double cd[2][2]={0};
    // double tpvx[8][8]={0};
    // double tpvy[8][8]={0};
    // double tpvu[8][8]={0};
    // double tpvv[8][8]={0};
    // struct disprm *dis=NULL;
    // double pv1[GAL_WCS_MAX_PVSIZE] = {0};
    // double pv2[GAL_WCS_MAX_PVSIZE] = {0};
    wcs=gal_wcs_read("test-pv.fits", "1", 0, 0, &a);

    struct disprm *dispre=NULL;
    dispre = malloc(sizeof(struct disprm));
    struct linprm lin;
    dispre->flag = -1;
    dispre=wcs->lin.dispre;
    lin=wcs->lin;
    lindist(1, &lin, dispre, -1);
    
    printf("%d\n", lin.dispre->naxis);
  return 0;
}
