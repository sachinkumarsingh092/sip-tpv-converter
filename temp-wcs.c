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
#include <wcslib/wcs.h>
#include <wcslib/wcslib.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/permutation.h>

static void
gal_wcs_calc_tpveq(double cd[2][2], double tpvu[8][8], double tpvv[8][8], \
                   char *infile, char *inhdu);
static double
gal_wcs_calcsip(size_t axis, size_t m, size_t n, \
                double tpvu[8][8], double tpvv[8][8]);



/***********************************/

#define GAL_WCS_MAX_PVSIZE    40
#define GAL_WCS_MAX_POLYORDER  8

#define max(a,b)                \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b)                \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })



/***********************************/







/* add filename etc. parameters later.*/
static void
gal_wcs_get_tpvparams(struct wcsprm *wcs,
                      double cd[2][2],
                      double *pv1,
                      double *pv2)
{
  size_t i, j, k,index=0;
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

  for(j=0; j < wcs->npvmax; j++)
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

        /* For ==> PVi_%dn", index, pv1[index]);
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
gal_wcs_get_sipparams(struct wcsprm *wcs,
                      double cd[2][2],
                      double a_coeff[5][5],
                      double b_coeff[5][5])
{
  const char *cp;
  double keyval=0;
  double *temp_cd=NULL;
  struct dpkey  *keyp=NULL;
  struct disprm *dispre=NULL;
  size_t i=0, j=0, m=0, n=0, k=0;

  dispre=wcs->lin.dispre;
  keyp=dispre->dp;

  temp_cd=gal_wcs_warp_matrix(wcs);


  for(i=0,k=0; i<2; i++)
    for(j=0; j<2; j++)
      {
        /* If a value is present store it, else store 0.*/
        if(temp_cd[k] != 0)   cd[i][j]=temp_cd[k++];
        else                  cd[i][j]=0;

        /* For a check:
        printf("CD%ld_%ld\t%.10E\n", i+1, j+1, cd[i][j]);
        */
      }


  for (i = 0; i < dispre->ndp; i++, keyp++)
    {
      if (keyp->j == 1)
        {
          if ((keyval = dpkeyd(keyp)) == 0.0) continue;

          cp = strchr(keyp->field, '.') + 1;
          if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
          cp += 8;

          sscanf(cp, "%ld_%ld", &m, &n);

          a_coeff[m][n]=dispre->dp[i].value.f;

          /*For a check.
          printf("A_%ld_%ld =\t %.10E\n", m, n, dispre->dp[i].value.f);
          */
        }
      else if (keyp->j == 2)
        {
          if ((keyval = dpkeyd(keyp)) == 0.0) continue;

          cp = strchr(keyp->field, '.') + 1;
          if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
          cp += 8;

          sscanf(cp, "%ld_%ld", &m, &n);

          b_coeff[m][n]=dispre->dp[i].value.f;

          /*For a check.
          printf("B_%ld_%ld =\t %.10E\n", m, n, dispre->dp[i].value.f);
          */
        }
    }
}





/*char *infile, inhdu in, naxis, npix to be
 transfered in main sip wrirting function from here.*/

static void
gal_wcs_get_sipcoeff(struct wcsprm *wcs,
                      size_t *a_order,
                      size_t *b_order,
                      double a_coeff[5][5],
                      double b_coeff[5][5],
                      char *infile,
                      char *inhdu)
{
  size_t m, n;
  double val=0;
  size_t i, j, k;
  double cd[2][2]={0};
  double tpvu[8][8]={0}, tpvv[8][8]={0};


  gal_wcs_calc_tpveq(cd, tpvu, tpvv, infile, inhdu);


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

  for(m=0; m<=4; m++)
    for(n=0; n<=4; n++)
      {
        /*For axis = 1*/
        val=gal_wcs_calcsip(1, m, n, tpvu, tpvv);
        if(val != 0)
          {
            a_coeff[m][n]=val;
            *a_order=max(*a_order, max(m,n));
          }
        else
          a_coeff[m][n]=0;

        /*For axis = 2*/
        val=gal_wcs_calcsip(2, m, n, tpvu, tpvv);
        if(val != 0)
          {
            b_coeff[m][n]=val;
            *b_order=max(*b_order, max(m,n));
          }
        else
          b_coeff[m][n]=0;
      }

  /*For a check.
  for(m=0;m<=4;m++)
    for(n=0;n<=4;n++)
      printf("A_%ld_%ld = %.10E \t B_%ld_%ld = %.10E\n",
              m, n, a_coeff[m][n], m, n, b_coeff[m][n]);
  */
}







/* Maybe tpv* matrix maybe not needed for this task.
   Remove it if redundant. */
static void
gal_wcs_intermidate_tpveq(double cd[2][2],
                    double *pv1,
                    double *pv2,
                    double k[5][5],
                    double l[5][5])
{
  size_t i, j;

  /* The indexes of the tpv* matrices are the exponents of x and y in
     a TPV distortion equation leaving out radial terms
     PV[1,3], PV[1,11], PV[1,23], PV[1,39] */

  /* Intermidate polynomials, k[i][j] and l[i][j].  */

  k[0][0]=pv1[0];
  l[0][0]=pv2[0];

  k[0][1]=cd[0][1]*pv1[1]+cd[1][1]*pv1[2];
  l[0][1]=cd[0][1]*pv2[2]+cd[1][1]*pv2[1];

  k[1][0]=cd[0][0]*pv1[1]+cd[1][0]*pv1[2];
  l[1][0]=cd[0][0]*pv2[2]+cd[1][0]*pv2[1];

  k[0][2]=cd[0][1]*cd[0][1]*pv1[4]+ \
          cd[0][1]*cd[1][1]*pv1[5]+ \
          cd[1][1]*cd[1][1]*pv1[6];
  l[0][2]=cd[0][1]*cd[0][1]*pv2[6]+ \
          cd[0][1]*cd[1][1]*pv2[5]+ \
          cd[1][1]*cd[1][1]*pv2[4];

  k[1][1]=2*cd[0][0]*cd[0][1]*pv1[4]+ \
            cd[0][0]*cd[1][1]*pv1[5]+ \
            cd[0][1]*cd[1][0]*pv1[5]+ \
          2*cd[1][0]*cd[1][1]*pv1[6];
  l[1][1]=2*cd[0][0]*cd[0][1]*pv2[6]+ \
            cd[0][0]*cd[1][1]*pv2[5]+ \
            cd[0][1]*cd[1][0]*pv2[5]+ \
          2*cd[1][0]*cd[1][1]*pv2[4];

  k[2][0]=cd[0][0]*cd[0][0]*pv1[4]+ \
          cd[0][0]*cd[1][0]*pv1[5]+ \
          cd[1][0]*cd[1][0]*pv1[6];
  l[2][0]=cd[0][0]*cd[0][0]*pv2[6]+ \
          cd[0][0]*cd[1][0]*pv2[5]+ \
          cd[1][0]*cd[1][0]*pv2[4];

  k[0][3]=cd[0][1]*cd[0][1]*cd[0][1]*pv1[7]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*pv1[8]+ \
          cd[0][1]*cd[1][1]*cd[1][1]*pv1[9]+ \
          cd[1][1]*cd[1][1]*cd[1][1]*pv1[10];
  l[0][3]=cd[0][1]*cd[0][1]*cd[0][1]*pv2[10]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*pv2[9]+ \
          cd[0][1]*cd[1][1]*cd[1][1]*pv2[8]+ \
          cd[1][1]*cd[1][1]*cd[1][1]*pv2[7];


  k[1][2]=3*cd[0][0]*cd[0][1]*cd[0][1]*pv1[7]+ \
          2*cd[0][0]*cd[0][1]*cd[1][1]*pv1[8]+ \
            cd[0][0]*cd[1][1]*cd[1][1]*pv1[9]+ \
            cd[0][1]*cd[0][1]*cd[1][0]*pv1[8]+ \
          2*cd[0][1]*cd[1][0]*cd[1][1]*pv1[9]+ \
          3*cd[1][0]*cd[1][1]*cd[1][1]*pv1[10];
  l[1][2]=3*cd[0][0]*cd[0][1]*cd[0][1]*pv2[10]+ \
          2*cd[0][0]*cd[0][1]*cd[1][1]*pv2[9]+ \
            cd[0][0]*cd[1][1]*cd[1][1]*pv2[8]+ \
            cd[0][1]*cd[0][1]*cd[1][0]*pv2[9]+ \
          2*cd[0][1]*cd[1][0]*cd[1][1]*pv2[8]+ \
          3*cd[1][0]*cd[1][1]*cd[1][1]*pv2[7];


  k[2][1]=3*cd[0][0]*cd[0][0]*cd[0][1]*pv1[7]+ \
          2*cd[0][0]*cd[0][1]*cd[1][0]*pv1[8]+ \
            cd[0][0]*cd[0][0]*cd[1][1]*pv1[8]+ \
            cd[0][1]*cd[1][0]*cd[1][0]*pv1[9]+ \
          2*cd[0][0]*cd[1][0]*cd[1][1]*pv1[9]+ \
          3*cd[1][0]*cd[1][0]*cd[1][1]*pv1[10];
  l[2][1]=3*cd[0][0]*cd[0][0]*cd[0][1]*pv2[10]+ \
          2*cd[0][0]*cd[0][1]*cd[1][0]*pv2[9]+ \
            cd[0][0]*cd[0][0]*cd[1][1]*pv2[9]+ \
            cd[0][1]*cd[1][0]*cd[1][0]*pv2[8]+ \
          2*cd[0][0]*cd[1][0]*cd[1][1]*pv2[8]+ \
          3*cd[1][0]*cd[1][0]*cd[1][1]*pv2[7];


  k[3][0]=cd[0][0]*cd[0][0]*cd[0][0]*pv1[7]+ \
          cd[0][0]*cd[0][0]*cd[1][0]*pv1[8]+ \
          cd[0][0]*cd[1][0]*cd[1][0]*pv1[9]+ \
          cd[1][0]*cd[1][0]*cd[1][0]*pv1[10];
  l[3][0]=cd[0][0]*cd[0][0]*cd[0][0]*pv2[10]+ \
          cd[0][0]*cd[0][0]*cd[1][0]*pv2[9]+ \
          cd[0][0]*cd[1][0]*cd[1][0]*pv2[8]+ \
          cd[1][0]*cd[1][0]*cd[1][0]*pv2[7];


  k[0][4]=cd[0][1]*cd[0][1]*cd[0][1]*cd[0][1]*pv1[12]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*cd[1][1]*pv1[13]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*cd[1][1]*pv1[14]+ \
          cd[0][1]*cd[1][1]*cd[1][1]*cd[1][1]*pv1[15]+ \
          cd[1][1]*cd[1][1]*cd[1][1]*cd[1][1]*pv1[16];
  l[0][4]=cd[0][1]*cd[0][1]*cd[0][1]*cd[0][1]*pv2[16]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*cd[1][1]*pv2[15]+ \
          cd[0][1]*cd[0][1]*cd[1][1]*cd[1][1]*pv2[14]+ \
          cd[0][1]*cd[1][1]*cd[1][1]*cd[1][1]*pv2[13]+ \
          cd[1][1]*cd[1][1]*cd[1][1]*cd[1][1]*pv2[12];


  k[1][3]=4*cd[0][0]*cd[0][1]*cd[0][1]*cd[0][1]*pv1[12]+ \
          3*cd[0][0]*cd[0][1]*cd[0][1]*cd[1][1]*pv1[13]+ \
          2*cd[0][0]*cd[0][1]*cd[1][1]*cd[1][1]*pv1[14]+ \
            cd[0][0]*cd[1][1]*cd[1][1]*cd[1][1]*pv1[15]+ \
            cd[0][1]*cd[0][1]*cd[0][1]*cd[1][0]*pv1[13]+ \
          2*cd[0][1]*cd[0][1]*cd[1][0]*cd[1][1]*pv1[14]+ \
          3*cd[0][1]*cd[1][0]*cd[1][1]*cd[1][1]*pv1[15]+ \
          4*cd[1][0]*cd[1][1]*cd[1][1]*cd[1][1]*pv1[16];
  l[1][3]=4*cd[0][0]*cd[0][1]*cd[0][1]*cd[0][1]*pv2[16]+ \
          3*cd[0][0]*cd[0][1]*cd[0][1]*cd[1][1]*pv2[15]+ \
          2*cd[0][0]*cd[0][1]*cd[1][1]*cd[1][1]*pv2[14]+ \
            cd[0][0]*cd[1][1]*cd[1][1]*cd[1][1]*pv2[13]+ \
            cd[0][1]*cd[0][1]*cd[0][1]*cd[1][0]*pv2[15]+ \
          2*cd[0][1]*cd[0][1]*cd[1][0]*cd[1][1]*pv2[14]+ \
          3*cd[0][1]*cd[1][0]*cd[1][1]*cd[1][1]*pv2[13]+ \
          4*cd[1][0]*cd[1][1]*cd[1][1]*cd[1][1]*pv2[12];


  k[2][2]=6*cd[0][0]*cd[0][0]*cd[0][1]*cd[0][1]*pv1[12]+ \
          3*cd[0][0]*cd[0][0]*cd[0][1]*cd[1][1]*pv1[13]+ \
            cd[0][0]*cd[0][0]*cd[1][1]*cd[1][1]*pv1[14]+ \
          3*cd[0][0]*cd[0][1]*cd[0][1]*cd[1][0]*pv1[13]+ \
          4*cd[0][0]*cd[0][1]*cd[1][0]*cd[1][1]*pv1[14]+ \
          3*cd[0][0]*cd[1][0]*cd[1][1]*cd[1][1]*pv1[15]+ \
            cd[0][1]*cd[0][1]*cd[1][0]*cd[1][0]*pv1[14]+ \
          3*cd[0][1]*cd[1][0]*cd[1][0]*cd[1][1]*pv1[15]+ \
          6*cd[1][0]*cd[1][0]*cd[1][1]*cd[1][1]*pv1[16];
  l[2][2]=6*cd[0][0]*cd[0][0]*cd[0][1]*cd[0][1]*pv2[16]+ \
          3*cd[0][0]*cd[0][0]*cd[0][1]*cd[1][1]*pv2[15]+ \
            cd[0][0]*cd[0][0]*cd[1][1]*cd[1][1]*pv2[14]+ \
          3*cd[0][0]*cd[0][1]*cd[0][1]*cd[1][0]*pv2[15]+ \
          4*cd[0][0]*cd[0][1]*cd[1][0]*cd[1][1]*pv2[14]+ \
          3*cd[0][0]*cd[1][0]*cd[1][1]*cd[1][1]*pv2[15]+ \
            cd[0][1]*cd[0][1]*cd[1][0]*cd[1][0]*pv2[14]+ \
          3*cd[0][1]*cd[1][0]*cd[1][0]*cd[1][1]*pv2[13]+ \
          6*cd[1][0]*cd[1][0]*cd[1][1]*cd[1][1]*pv2[12];

  k[3][1]=4*cd[0][0]*cd[0][0]*cd[0][0]*cd[0][1]*pv1[12]+ \
            cd[0][0]*cd[0][0]*cd[0][0]*cd[1][1]*pv1[13]+ \
          3*cd[0][0]*cd[0][0]*cd[0][1]*cd[1][0]*pv1[13]+ \
          2*cd[0][0]*cd[0][0]*cd[1][0]*cd[1][1]*pv1[14]+ \
          2*cd[0][0]*cd[0][1]*cd[1][0]*cd[1][0]*pv1[14]+ \
          3*cd[0][0]*cd[1][0]*cd[1][0]*cd[1][1]*pv1[15]+ \
            cd[0][1]*cd[1][0]*cd[1][0]*cd[1][0]*pv1[15]+ \
          4*cd[1][0]*cd[1][0]*cd[1][0]*cd[1][1]*pv1[16];
  l[3][1]=4*cd[0][0]*cd[0][0]*cd[0][0]*cd[0][1]*pv2[16]+ \
            cd[0][0]*cd[0][0]*cd[0][0]*cd[1][1]*pv2[15]+ \
          3*cd[0][0]*cd[0][0]*cd[0][1]*cd[1][0]*pv2[15]+ \
          2*cd[0][0]*cd[0][0]*cd[1][0]*cd[1][1]*pv2[14]+ \
          2*cd[0][0]*cd[0][1]*cd[1][0]*cd[1][0]*pv2[14]+ \
          3*cd[0][0]*cd[1][0]*cd[1][0]*cd[1][1]*pv2[13]+ \
            cd[0][1]*cd[1][0]*cd[1][0]*cd[1][0]*pv2[13]+ \
          4*cd[1][0]*cd[1][0]*cd[1][0]*cd[1][1]*pv2[12];


  k[4][0]=cd[0][0]*cd[0][0]*cd[0][0]*cd[0][0]*pv1[12]+ \
          cd[0][0]*cd[0][0]*cd[0][0]*cd[1][0]*pv1[13]+ \
          cd[0][0]*cd[0][0]*cd[1][0]*cd[1][0]*pv1[14]+ \
          cd[0][0]*cd[1][0]*cd[1][0]*cd[1][0]*pv1[15]+ \
          cd[1][0]*cd[1][0]*cd[1][0]*cd[1][0]*pv1[16];
  l[4][0]=cd[0][0]*cd[0][0]*cd[0][0]*cd[0][0]*pv2[16]+ \
          cd[0][0]*cd[0][0]*cd[0][0]*cd[1][0]*pv2[15]+ \
          cd[0][0]*cd[0][0]*cd[1][0]*cd[1][0]*pv2[14]+ \
          cd[0][0]*cd[1][0]*cd[1][0]*cd[1][0]*pv2[13]+ \
          cd[1][0]*cd[1][0]*cd[1][0]*cd[1][0]*pv2[12];



  /* For a check:
  for(i=0; i<=4; i++)
    for(j=0;j<=4;j++)
      {
        printf("k%ld_%ld \t %.10E\n",   i, j, k[i][j]);
        printf("l%ld_%ld \t %.10E\n\n", i, j, l[i][j]);
      }
  */

}




static void
gal_wcs_intermidate_sipeq(double cd[2][2],
                          double cd_inv[2][2],
                          double a_coeff[5][5],
                          double b_coeff[5][5],
                          double *pv1,
                          double *pv2)
{
  size_t i, j;

  /* Intemidiate equations which give the value of PV.
     First index is the axis and second index is the coefficient
     number in the original tpv* equation. */

    // pvi_0  
  pv1[0]=a_coeff[0][0]*cd[0][0] + b_coeff[0][0]*cd[0][1];
  pv2[0]=a_coeff[0][0]*cd[1][0] + b_coeff[0][0]*cd[1][1];


    // pvi_1  
  pv1[1]=a_coeff[0][1]*cd[0][0]*cd_inv[1][0]+\
        a_coeff[1][0]*cd[0][0]*cd_inv[0][0]+\
        b_coeff[0][1]*cd[0][1]*cd_inv[1][0]+\
        b_coeff[1][0]*cd[0][1]*cd_inv[0][0]+\
        cd[0][0]*cd_inv[0][0]+\
        cd[0][1]*cd_inv[1][0];
  pv2[1]=a_coeff[0][1]*cd[1][0]*cd_inv[1][1]+\
        a_coeff[1][0]*cd[1][0]*cd_inv[0][1]+\
        b_coeff[0][1]*cd[1][1]*cd_inv[1][1]+\
        b_coeff[1][0]*cd[1][1]*cd_inv[0][1]+\
        cd[1][0]*cd_inv[0][1]+\
        cd[1][1]*cd_inv[1][1];



    //  pvi_2  
  pv1[2]=a_coeff[0][1]*cd[0][0]*cd_inv[1][1]+\
        a_coeff[1][0]*cd[0][0]*cd_inv[0][1]+\
        b_coeff[0][1]*cd[0][1]*cd_inv[1][1]+\
        b_coeff[1][0]*cd[0][1]*cd_inv[0][1]+\
        cd[0][0]*cd_inv[0][1]+\
        cd[0][1]*cd_inv[1][1];
  pv2[2]=a_coeff[0][1]*cd[1][0]*cd_inv[1][0]+\
        a_coeff[1][0]*cd[1][0]*cd_inv[0][0]+\
        b_coeff[0][1]*cd[1][1]*cd_inv[1][0]+\
        b_coeff[1][0]*cd[1][1]*cd_inv[0][0]+\
        cd[1][0]*cd_inv[0][0]+\
        cd[1][1]*cd_inv[1][0];



    // pvi_4  
  pv1[4]=a_coeff[0][2]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]+\
        a_coeff[1][1]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]+\
        a_coeff[2][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]+\
        b_coeff[0][2]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]+\
        b_coeff[1][1]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]+\
        b_coeff[2][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0];
  pv2[4]=a_coeff[0][2]*cd[1][0]*cd_inv[1][1]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[1][0]*cd_inv[0][1]*cd_inv[1][1]+\
        a_coeff[2][0]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]+\
        b_coeff[0][2]*cd[1][1]*cd_inv[1][1]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[1][1]*cd_inv[0][1]*cd_inv[1][1]+\
        b_coeff[2][0]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1];



    //  pvi_5  
  pv1[5]=2*a_coeff[0][2]*cd[0][0]*cd_inv[1][0]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[0][0]*cd_inv[0][0]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]+\
        2*a_coeff[2][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]+\
        2*b_coeff[0][2]*cd[0][1]*cd_inv[1][0]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[0][1]*cd_inv[0][0]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]+\
        2*b_coeff[2][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1];
  pv2[5]=2*a_coeff[0][2]*cd[1][0]*cd_inv[1][0]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[1][0]*cd_inv[0][0]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]+\
        2*a_coeff[2][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]+\
        2*b_coeff[0][2]*cd[1][1]*cd_inv[1][0]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[1][1]*cd_inv[0][0]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]+\
        2*b_coeff[2][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1];



    //  pvi_6  
  pv1[6]=a_coeff[0][2]*cd[0][0]*cd_inv[1][1]*cd_inv[1][1]+\
        a_coeff[1][1]*cd[0][0]*cd_inv[0][1]*cd_inv[1][1]+\
        a_coeff[2][0]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]+\
        b_coeff[0][2]*cd[0][1]*cd_inv[1][1]*cd_inv[1][1]+\
        b_coeff[1][1]*cd[0][1]*cd_inv[0][1]*cd_inv[1][1]+\
        b_coeff[2][0]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1];
  pv2[6]=a_coeff[0][2]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]+\
        a_coeff[1][1]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]+\
        a_coeff[2][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]+\
        b_coeff[0][2]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]+\
        b_coeff[1][1]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]+\
        b_coeff[2][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0];



    //  pvi_7  
  pv1[7]=a_coeff[0][3]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]+\
        a_coeff[1][2]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]+\
        a_coeff[2][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]+\
        a_coeff[3][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]+\
        b_coeff[0][3]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]+\
        b_coeff[1][2]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]+\
        b_coeff[2][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]+\
        b_coeff[3][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0];
  pv2[7]=a_coeff[0][3]*cd[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]+\
        a_coeff[1][2]*cd[1][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]+\
        a_coeff[2][1]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]+\
        a_coeff[3][0]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]+\
        b_coeff[0][3]*cd[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]+\
        b_coeff[1][2]*cd[1][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]+\
        b_coeff[2][1]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]+\
        b_coeff[3][0]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];


    //  pvi_8  
  pv1[8]=3*a_coeff[0][3]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1]+\
         2*a_coeff[1][2]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1]+\
         a_coeff[1][2]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]+\
         a_coeff[2][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
         2*a_coeff[2][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
         3*a_coeff[3][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1] +\
         3*b_coeff[0][3]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
         2*b_coeff[1][2]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
         b_coeff[1][2]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
         b_coeff[2][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
         2*b_coeff[2][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
         3*b_coeff[3][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1];
   pv2[8]=3*a_coeff[0][3]*cd[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][2]*cd[1][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] + \
          2*a_coeff[1][2]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          2*a_coeff[2][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[2][1]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          3*a_coeff[3][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1] +\
          3*b_coeff[0][3]*cd[1][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][2]*cd[1][1]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[1][2]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          2*b_coeff[2][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[2][1]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          3*b_coeff[3][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1];




    //  pvi_9  
  pv1[9]=3*a_coeff[0][3]*cd[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][2]*cd[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*a_coeff[1][2]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          2*a_coeff[2][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[2][1]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          3*a_coeff[3][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1] +\
          3*b_coeff[0][3]*cd[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][2]*cd[0][1]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[1][2]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          2*b_coeff[2][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[2][1]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          3*b_coeff[3][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1];

  pv2[9]=3*a_coeff[0][3]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*a_coeff[1][2]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[1][2]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[2][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          2*a_coeff[2][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          3*a_coeff[3][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1] +\
          3*b_coeff[0][3]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*b_coeff[1][2]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[1][2]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[2][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          2*b_coeff[2][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          3*b_coeff[3][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1];



    // pvi_10  
  pv1[10]=a_coeff[0][3]*cd[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][2]*cd[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[2][1]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[3][0]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1] +\
          b_coeff[0][3]*cd[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][2]*cd[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[2][1]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[3][0]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];


  pv2[10]=a_coeff[0][3]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[1][2]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[2][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          a_coeff[3][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0] +\
          b_coeff[0][3]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[1][2]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[2][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          b_coeff[3][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0];





    //  pvi_12
  pv1[12]=a_coeff[0][4]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[1][3]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          a_coeff[4][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0] +\
          b_coeff[0][4]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[1][3]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          b_coeff[4][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0];

  pv2[12]=a_coeff[0][4]*cd[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[1][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[3][1]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[4][0]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1] +\
          b_coeff[0][4]*cd[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[1][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[3][1]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[4][0]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];




    //  pvi_13
  pv1[13]=4*a_coeff[0][4]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          2*a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          4*a_coeff[4][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1] +\
          4*b_coeff[0][4]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          2*b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          4*b_coeff[4][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1];


  pv2[13]=4*a_coeff[0][4]*cd[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[1][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[3][1]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          4*a_coeff[4][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1] +\
          4*b_coeff[0][4]*cd[1][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[1][1]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[3][1]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          4*b_coeff[4][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];



    //  pvi_14>
  pv1[14]=6*a_coeff[0][4]*cd[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          4*a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          3*a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          6*a_coeff[4][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1] +\
          6*b_coeff[0][4]*cd[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[0][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          4*b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          3*b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          6*b_coeff[4][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1];

  pv2[14]=6*a_coeff[0][4]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          4*a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[1][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          3*a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          6*a_coeff[4][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1] +\
          6*b_coeff[0][4]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1] +\
          4*b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[1][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          3*b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          6*b_coeff[4][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1];



    //  pvi_15> 
  pv1[15]=4*a_coeff[0][4]*cd[0][0]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[0][0]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[3][1]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          4*a_coeff[4][0]*cd[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1] +\
          4*b_coeff[0][4]*cd[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[0][1]*cd_inv[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[3][1]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][0] +\
          4*b_coeff[4][0]*cd[0][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];

  pv2[15]=4*a_coeff[0][4]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          3*a_coeff[1][3]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[1][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          2*a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          3*a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          4*a_coeff[4][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1] +\
          4*b_coeff[0][4]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          3*b_coeff[1][3]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[1][1]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          2*b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][1] +\
          2*b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][1] +\
          3*b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1]*cd_inv[1][0] +\
          4*b_coeff[4][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][1];



    //  pvi_16> 
  pv1[16]=a_coeff[0][4]*cd[0][0]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[1][3]*cd[0][0]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[2][2]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          a_coeff[3][1]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          a_coeff[4][0]*cd[0][0]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1] +\
          b_coeff[0][4]*cd[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[1][3]*cd[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[2][2]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1]*cd_inv[1][1] +\
          b_coeff[3][1]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[1][1] +\
          b_coeff[4][0]*cd[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1]*cd_inv[0][1];

  pv2[16]=a_coeff[0][4]*cd[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[1][3]*cd[1][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[2][2]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          a_coeff[3][1]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          a_coeff[4][0]*cd[1][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0] +\
          b_coeff[0][4]*cd[1][1]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[1][3]*cd[1][1]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[2][2]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0]*cd_inv[1][0] +\
          b_coeff[3][1]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[1][0] +\
          b_coeff[4][0]*cd[1][1]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0]*cd_inv[0][0];




    

  // /* For a check:
    for(j=0; j<=16 ;j++)
      {
        if (j == 3 || j == 11) continue;
        printf("pv1_%ld \t %.12E\n", j, pv1[j]);
        printf("pv2_%ld \t %.12E\n", j, pv2[j]);
      }
  // */


}




static void
gal_wcs_calc_tpveq(double cd[2][2],
                   double tpvu[8][8],
                   double tpvv[8][8],
                   char *infile,
                   char *inhdu)
{
  /* tpvu and tpvv are u-v distortion equations in TPV convention. */
  int a;
  size_t i,j;
  double determinant=0;
  struct wcsprm *wcs=NULL;
  double cd_inv[2][2]={0};
  double pv1[GAL_WCS_MAX_PVSIZE] = {0};
  double pv2[GAL_WCS_MAX_PVSIZE] = {0};
  double k[5][5]={0}, l[5][5]={0};
  // double tpv1[8][8]={0}, tpv2[8][8]={0};


  wcs=gal_wcs_read(infile, inhdu, 0, 0, &a);

  gal_wcs_get_tpvparams(wcs, cd, pv1, pv2);

  gal_wcs_intermidate_tpveq(cd, pv1, pv2, k, l);


  /* We will find matrix tpvu and tpvv by finding inverse of
     CD matrix and multiplying with tpv* matrix.
     For inverse of a 2x2 matrix we use the below trasformations:
              inverse(|a  b|) =  1 *|d  -b|
                      |c  d|    |A| |-c  a|
      where |A| is the determinant of the matrix which is calculated by:
                          |A| = a*d-b*c.
    */


  determinant = cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0];

  /* Inverse matrix */
  cd_inv[0][0]=cd[1][1]/determinant;      /* a */
  cd_inv[0][1]=(-1*cd[0][1])/determinant; /* b */
  cd_inv[1][0]=(-1*cd[1][0])/determinant; /* c */
  cd_inv[1][1]=cd[0][0]/determinant;      /* d */

  /*For a check.
  printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", cd_inv[0][0], cd_inv[0][1], \
                                             cd_inv[1][0], cd_inv[1][1]);
  printf("%.10lf\n", determinant);
  */

  /* For matrix tpvv and tpvu, we have to use the following
     matrix equation:

                  |tpvu| = cd_inv*|k[i][j]|
                  |tpvv|          |l[i][j]|
    though intermidate polynomial equations have to be calculated prior
    to multiplycation with cd_inv.
      */

  for(i=0; i<=4; i++)
    for(j=0; j<=4; j++){
      tpvu[i][j]=cd_inv[0][0]*k[i][j]+cd_inv[0][1]*l[i][j];
      tpvv[i][j]=cd_inv[1][0]*k[i][j]+cd_inv[1][1]*l[i][j];

      /*For a check:
      printf("%.10E, %.10E\n", tpvu[i][j], tpvv[i][j]);
      */
    }

}





static void
gal_wcs_calc_sipeq(double cd[2][2],
                   double *pv1,
                   double *pv2,
                   char *infile,
                   char *inhdu)
{
  /* tpvu and tpvv are u-v distortion equations in TPV convention. */
  int a;
  size_t i,j;
  double determinant=0;
  struct wcsprm *wcs=NULL;
  double cd_inv[2][2]={0};
  double a_coeff[5][5]={0}, b_coeff[5][5]={0};

  wcs=gal_wcs_read(infile, inhdu, 0, 0, &a);


  /* We will find matrix tpvu and tpvv by finding inverse of
     CD matrix and multiplying with tpv* matrix.
     For inverse of a 2x2 matrix we use the below trasformations:
              inverse(|a  b|) =  1 *|d  -b|
                      |c  d|    |A| |-c  a|
      where |A| is the determinant of the matrix which is calculated by:
                          |A| = a*d-b*c.
    */

  gal_wcs_get_sipparams(wcs, cd, a_coeff, b_coeff);
 
  determinant = cd[0][0]*cd[1][1] - cd[0][1]*cd[1][0];

  /* Inverse matrix */
  cd_inv[0][0]=cd[1][1]/determinant;      /* a */
  cd_inv[0][1]=(-1*cd[0][1])/determinant; /* b */
  cd_inv[1][0]=(-1*cd[1][0])/determinant; /* c */
  cd_inv[1][1]=cd[0][0]/determinant;      /* d */


  gal_wcs_intermidate_sipeq( cd, cd_inv, a_coeff, b_coeff, pv1, pv2);


  /*For a check.
  printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", cd_inv[0][0], cd_inv[0][1], \
                                             cd_inv[1][0], cd_inv[1][1]);
  printf("%.10lf\n", determinant);
  */


}



/* Calculate the SIP coefficients from CD matrix
    parameters and PV coefficients. */
static double
gal_wcs_calcsip(size_t axis,
                size_t m,
                size_t n,
                double tpvu[8][8],
                double tpvv[8][8])
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




/* Just combine fitreverse and evalpoly. */
static void
gal_wcs_fitreverse(double *u,
                 double *v,
                 size_t a_order,
                 size_t b_order,
                 size_t naxis1,
                 size_t naxis2,
                 double a_coeff[5][5],
                 double b_coeff[5][5],
                 double ap_coeff[5][5],
                 double bp_coeff[5][5])
{
  double chisq_ap, chisq_bp;
  double *udiff=NULL, *vdiff=NULL;
  double *uprime=NULL, *vprime=NULL;
  double **udict=NULL, **vdict=NULL;
  size_t tsize=(naxis1/4)*(naxis2/4);
  double **updict=NULL, **vpdict=NULL;
  gsl_vector *y_ap, *y_bp, *c_ap, *c_bp;
  size_t ap_order=a_order, bp_order=b_order;
  gsl_matrix *X_ap, *X_bp, *cov_ap, *cov_bp;
  size_t i=0, j=0, k=0, p_ap=0, p_bp=0, ij=0;


  /*Allocate updict and vpdict.*/
  updict=malloc((ap_order+1)*sizeof(*updict));
  vpdict=malloc((bp_order+1)*sizeof(*vpdict));

  for(i=0; i<=max(ap_order, bp_order); i++)
    {
      updict[i]=malloc(tsize*sizeof(updict));
      vpdict[i]=malloc(tsize*sizeof(vpdict));
    }


  /*allocating uprime and vprime.*/
  uprime=malloc(tsize*sizeof(*uprime));
  vprime=malloc(tsize*sizeof(*vprime));


  /*allocating uprime and vprime.*/
  udiff=malloc(tsize*sizeof(*udiff));
  vdiff=malloc(tsize*sizeof(*vdiff));


  /*allocating uprime and vprime.*/
  udict=malloc((a_order+1)*sizeof(*udict));
  vdict=malloc((b_order+1)*sizeof(*vdict));

  for(i=0; i<=max(a_order, b_order); i++)
    {
      udict[i]=malloc(tsize*sizeof(udict));
      vdict[i]=malloc(tsize*sizeof(vdict));
    }


  /* Initialize dicts. */
  for(i=0; i<tsize; i++)
    {
      udict[0][i]=1;
      vdict[0][i]=1;
    }


  /* Fill the values from the in the dicts. The rows of the
      dicts act as a key to achieve a key-value functionality. */
  for(i=1; i<=a_order; i++)
    for(j=0; j<tsize; j++)
      udict[i][j]=udict[i-1][j]*u[j];

  for(i=1; i<=b_order; i++)
    for(j=0; j<tsize; j++)
      vdict[i][j]=vdict[i-1][j]*v[j];

  /* Initialising uprime and vprime. */
  for(i=0; i<tsize; i++)
    {
      uprime[i]=u[i];
      vprime[i]=v[i];
    }


  /* Populating forward coefficients on a grid. */
  for(i=0; i<=a_order; i++)
    for(j=0; j<=a_order-i; j++)
      {
        for(k=0; k<tsize; k++)
          uprime[k]+=a_coeff[i][j]*udict[i][k]*vdict[j][k];

        /* The number of parameters for AP_* coefficients. */
        p_ap++;
      }


  for(i=0; i<=b_order; i++)
    for(j=0; j<=b_order-i; j++)
        {
          for(k=0; k<tsize; k++)
            vprime[k]+=b_coeff[i][j]*udict[i][k]*vdict[j][k];

          /* The number of parameters for BP_* coefficients. */
          p_bp++;
        }

  // printf("p=%ld\n", p_ap);
  /*For a check.
  for(i=0; i<=a_order; i++)
    for(j=0; j<tsize; j++){
      printf("udict[%ld][%ld] = %.8lf\n", i, j, uprime[j]);
      printf("u%ld = %.10E\n", i, u[i]);
  }
  */


  /***********evalpoly ends***************/


  /* Initialize dicts. */
  for(i=0; i<tsize; i++)
    {
      updict[0][i]=1;
      vpdict[0][i]=1;
    }


  /* Fill the values from the in the dicts. The rows of the
      dicts act as a key to achieve a key-value functionality. */
  for(i=1; i<=ap_order; i++)
    for(j=0; j<tsize; j++)
      updict[i][j]=updict[i-1][j]*uprime[j];

  for(i=1; i<=bp_order; i++)
    for(j=0; j<tsize; j++)
      vpdict[i][j]=vpdict[i-1][j]*vprime[j];


  /* Allocate memory for Multi-parameter Linear Regressions. */
  X_ap = gsl_matrix_alloc (tsize, p_ap);
  X_bp = gsl_matrix_alloc (tsize, p_bp);

  y_ap = gsl_vector_alloc (tsize);
  y_bp = gsl_vector_alloc (tsize);

  c_ap = gsl_vector_alloc (p_ap);
  c_bp = gsl_vector_alloc (p_bp);

  cov_ap = gsl_matrix_alloc (p_ap, p_ap);
  cov_bp = gsl_matrix_alloc (p_bp, p_bp);




  for(i=0; i<tsize; i++)
    {
      udiff[i]=u[i]-uprime[i];
      vdiff[i]=v[i]-vprime[i];

      gsl_vector_set (y_ap, i, udiff[i]);
      gsl_vector_set (y_bp, i, vdiff[i]);

    }


  /* Filling up he X_ap matrix in column-wise order. */
  for(i=0; i<=ap_order; i++)
    for(j=0; j<=ap_order-i; j++)
      {
        for(k=0; k<tsize; k++)
          {
            gsl_matrix_set (X_ap, k, ij, updict[i][k]*vpdict[j][k]);

            /*For a check.
            printf("x_ap[%ld] = %.8lf\n", ij, updict[i][k]*vpdict[j][k]);
            */
          }
        ij++;
      }


  ij=0;
  /* Filling up he X_bp matrix in column-wise order. */
  for(i=0; i<=bp_order; i++)
    for(j=0; j<=bp_order-i; j++)
      {
        for(k=0; k<tsize; k++)
          {
            gsl_matrix_set (X_bp, k, ij, updict[i][k]*vpdict[j][k]);

            /*For a check.
            printf("x_bp[%ld] = %.8lf\n", ij, updict[i][k]*vpdict[j][k]);
            */
          }
        ij++;
      }



  gsl_multifit_linear_workspace * work_ap =
        gsl_multifit_linear_alloc (tsize, p_ap);
  gsl_multifit_linear_workspace * work_bp =
        gsl_multifit_linear_alloc (tsize, p_bp);


  gsl_multifit_linear (X_ap, y_ap, c_ap, cov_ap,
                          &chisq_ap, work_ap);
  gsl_multifit_linear (X_bp, y_bp, c_bp, cov_bp,
                          &chisq_bp, work_bp);



  p_ap=0;
  for(i=0; i<=ap_order; i++)
    {
      for(j=0; j<=ap_order-i; j++)
        {
          ap_coeff[i][j]=gsl_vector_get(c_ap, p_ap);

          /*For a check.
          printf("AP_%ld_%ld = %.8E\n", i, j, ap_coeff[i][j]);
          */

          p_ap++;
        }
    }

  p_bp=0;
  for(i=0; i<=bp_order; i++)
    {
      for(j=0; j<=bp_order-i; j++)
        {
          bp_coeff[i][j]=gsl_vector_get(c_bp, p_bp);

          /*For a check.
          printf("BP_%ld_%ld = %.8E\n", i, j, bp_coeff[i][j]);
          */

          p_bp++;
        }
    }



  /*For a check.
  for(i=0; i<p_ap; i++)
    for(j=0; j<tsize; j++){
      printf("X[%ld][%ld] = %.8lf\n", i, j, gsl_matrix_get(X_ap, j, i));
  }
  */

  /* Free memory and return.*/
  gsl_multifit_linear_free (work_bp);
  gsl_multifit_linear_free (work_ap);

  gsl_matrix_free (cov_bp);
  gsl_matrix_free (cov_ap);

  gsl_vector_free (c_bp);
  gsl_vector_free (c_ap);

  gsl_vector_free (y_bp);
  gsl_vector_free (y_ap);

  gsl_matrix_free (X_bp);
  gsl_matrix_free (X_ap);


  for(i=0; i<a_order; i++)
    {
      free(vdict[i]);
      free(udict[i]);
    }
  free(vdict);
  free(udict);

  free(vdiff);
  free(udiff);

  free(vprime);
  free(uprime);


  for(i=0; i<ap_order; i++)
    {
      free(vpdict[i]);
      free(updict[i]);
    }
  free(vpdict);
  free(updict);

}



/* Add wcs to read crpix and naxisn here
    instead of read_sipparams.*/
void
gal_wcs_get_revkeyvalues(struct wcsprm *wcs,
                        gal_data_t *in,
                        double ap_coeff[5][5],
                        double bp_coeff[5][5],
                        char *infile,
                        char *inhdu)
{
  size_t tsize=0;
  size_t i=0, j=0, k=0;
  size_t naxis1, naxis2;
  double crpix1, crpix2;
  double *u=NULL, *v=NULL;
  size_t a_order=0, b_order=0;
  double a_coeff[5][5]={0}, b_coeff[5][5]={0};

  in=gal_fits_img_read(infile, inhdu, -1, 1);

  naxis2=in->dsize[0];
  naxis1=in->dsize[1];

  crpix1=wcs->crpix[0];
  crpix2=wcs->crpix[1];

  tsize=(naxis1/4)*(naxis2/4);


  /* Allocate the size of u,v arrays. */
  u=malloc(tsize*sizeof(*u));
  v=malloc(tsize*sizeof(*v));

  for(i=0; i<naxis2; i+=4)
    for(j=0; j<naxis1; j+=4)
      u[k++]=j-crpix1;

  k=0;
  for(i=0; i<naxis2; i+=4)
    for(j=0; j<naxis1; j+=4)
      v[k++]=i-crpix2;

  /*For a check.
  for(i=0; i<tsize; i++)
    printf("u%ld = %.10E\n", i, u[i]);
  */

  gal_wcs_get_sipcoeff(wcs, &a_order, &b_order, a_coeff, b_coeff, infile, inhdu);

  gal_wcs_fitreverse(u, v, a_order, b_order, naxis1, naxis2,
                     a_coeff, b_coeff, ap_coeff, bp_coeff);


  free(v);
  free(u);

}







char *
gal_wcs_add_sipkeywords(struct wcsprm *wcs,
                        gal_data_t *in,
                        double tpvu[8][8],
                        double tpvv[8][8],
                        double cd[2][2],
                        char *infile,
                        char *inhdu,
                        int add_reverse,
                        int *nkeys)
{
  double val=0;
  uint8_t i, j, k=0;
  int size = wcs->naxis;
  size_t a_order=0, b_order=0;
  size_t ap_order=0, bp_order=0;
  size_t m, n, num=0, numkey=100;
  double ap_coeff[5][5]={0}, bp_coeff[5][5]={0};
  char *fullheader, fmt[50], sipkey[8], keyaxis[9], pcaxis[10];

  *nkeys = 0;

  /* The format for each card. */
  sprintf(fmt, "%%-8s= %%20.12E%%50s");

  gal_wcs_calc_tpveq(cd, tpvu, tpvv, infile, inhdu);


  /* Allcate memory for cards. */
  fullheader = malloc(numkey*80);

  /* Add other necessary cards. */
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20d%50s", "WCSAXES", wcs->naxis, "");

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CRPIX%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.8lf%50s", keyaxis, wcs->crpix[i-1], "");
    }

  for(i=1; i<=size; i++)
    for(j=1; j<=size; j++)
      {
        sprintf(pcaxis, "PC%d_%d", i, j);
        sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", pcaxis, wcs->pc[k++], "");
      }

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CDELT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", keyaxis, wcs->cdelt[i-1], "");
    }

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CUNIT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", keyaxis, wcs->cunit[i-1]);
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE1", "'RA---TAN-SIP'");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE2", "'DEC--TAN-SIP'");

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CRVAL%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.10lf%50s", keyaxis, wcs->crval[i-1], "");
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LONPOLE", wcs->lonpole, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LATPOLE", wcs->latpole, "");

  for(i=1; i<=size; i++)
    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "MJDREFI", wcs->mjdref[i-1], "");

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "RADESYS", wcs->radesys);

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "EQUINOX", wcs->equinox, "");


  for(m=0; m<=4; m++)
    for(n=0; n<=4; n++)
      {
        /*For axis = 1*/
        val=gal_wcs_calcsip(1, m, n, tpvu, tpvv);
        if(val != 0)
          {
            /* Make keywords */
            sprintf(sipkey, "A_%ld_%ld", m, n);
            sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
            a_order=max(a_order, max(m,n));
          }

        /*For axis = 2*/
        val=gal_wcs_calcsip(2, m, n, tpvu, tpvv);
        if(val != 0)
          {
            /* Make keywords */
            sprintf(sipkey, "B_%ld_%ld", m, n);
            sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
            b_order=max(b_order, max(m,n));
          }

      }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "A_ORDER", a_order, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "B_ORDER", b_order, "");


  if( add_reverse )
    {
      ap_order=a_order;
      bp_order=b_order;

      gal_wcs_get_revkeyvalues(wcs, in, ap_coeff, bp_coeff, infile, inhdu);

      for(m=0; m<=ap_order; m++)
        for(n=0; n<=ap_order-m; n++)
          {
            /*For axis = 1*/
            val=ap_coeff[m][n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(sipkey, "AP_%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
              }

            /*For axis = 2*/
            val=bp_coeff[m][n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(sipkey, "BP_%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, sipkey, val, "");
              }
          }

    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "AP_ORDER", ap_order, "");
    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20ld%50s", "BP_ORDER", bp_order, "");

    }



  *nkeys = num;

  /*For a check.
  printf("%s\n", fullheader);
  */

  return fullheader;

}





char *
gal_wcs_add_pvkeywords(struct wcsprm *wcs,
                        gal_data_t *in,
                        double *pv1,
                        double *pv2,
                        double cd[2][2],
                        char *infile,
                        char *inhdu,
                        int *nkeys)
{
  double val=0;
  uint8_t i, j, k=0;
  int size = wcs->naxis;
  size_t m, n, num=0, numkey=100;
  char *fullheader, fmt[50], pvkey[8], keyaxis[9], pcaxis[10];

  *nkeys = 0;

  /* The format for each card. */
  sprintf(fmt, "%%-8s= %%20.12E%%50s");

  gal_wcs_calc_sipeq(cd, pv1, pv2, infile, inhdu);


  /* Allcate memory for cards. */
  fullheader = malloc(numkey*80);

  /* Add other necessary cards. */
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20d%50s", "WCSAXES", wcs->naxis, "");

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CRPIX%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.8lf%50s", keyaxis, wcs->crpix[i-1], "");
    }

  for(i=1; i<=size; i++)
    for(j=1; j<=size; j++)
      {
        sprintf(pcaxis, "PC%d_%d", i, j);
        sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", pcaxis, wcs->pc[k++], "");
      }

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CDELT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", keyaxis, wcs->cdelt[i-1], "");
    }

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CUNIT%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", keyaxis, wcs->cunit[i-1]);
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE1", "'RA---TPV'");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "CTYPE2", "'DEC--TPV'");

  for(i=1; i<=size; i++)
    {
      sprintf(keyaxis, "CRVAL%d", i);
      sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.10lf%50s", keyaxis, wcs->crval[i-1], "");
    }

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LONPOLE", wcs->lonpole, "");
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.17lf%50s", "LATPOLE", wcs->latpole, "");

  for(i=1; i<=size; i++)
    sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "MJDREFI", wcs->mjdref[i-1], "");

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %-70s", "RADESYS", wcs->radesys);

  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20.1lf%50s", "EQUINOX", wcs->equinox, "");

  for(m=1; m<=2; m++)
    for(n=0; n<=16; n++)
      {
        /*For axis = 1*/
        if(m == 1)
          {
            val=pv1[n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(pvkey, "PV%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, pvkey, val, "");
              }
          }

        /*For axis = 2*/
        if(m == 2)
          {
            val=pv2[n];
            if(val != 0)
              {
                /* Make keywords */
                sprintf(pvkey, "PV%ld_%ld", m, n);
                sprintf(fullheader+(FLEN_CARD-1)*num++, fmt, pvkey, val, "");
              }
          }

      }


  *nkeys = num;

  /*For a check.
  printf("%s\n", fullheader);
  */

  return fullheader;

}





struct wcsprm *
gal_wcs_tpv2sip(struct wcsprm *inwcs,
                gal_data_t *in,
                char *infile,
                char *inhdu)
{
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */
  size_t i, fulllen;
  int nwcs, sumcheck;
  int nkeys=0, status=0;
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  struct wcsprm *outwcs=NULL;
  double cd[2][2]={0};
  double tpvu[8][8] ={0}, tpvv[8][8]={0};

  char *fullheader=gal_wcs_add_sipkeywords(inwcs, in, tpvu, tpvv, cd, infile, inhdu, 1, &nkeys);


  // printf("%s\n", fullheader);
  /* WCSlib function to parse the FITS headers. */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, &nwcs, &outwcs);
  if( outwcs == NULL )
  {

    fprintf(stderr, "\n##################\n"
            "WCSLIB Warning: wcspih ERROR %d: %s.\n"
            "##################\n",
            status, wcs_errmsg[status]);
    outwcs=NULL; nwcs=0;
  }

  /* Set the internal structure: */
  if(outwcs)
    {
      /* It may happen that the WCS-related keyword values are stored as
         strings (they have single-quotes around them). In this case,
         WCSLIB will read the CRPIX and CRVAL values as zero. When this
         happens do a small check and abort, while informing the user about
         the problem. */
      sumcheck=0;
      for(i=0;i<outwcs->naxis;++i)
        {sumcheck += (outwcs->crval[i]==0.0f) + (outwcs->crpix[i]==0.0f);}
      if(sumcheck==outwcs->naxis*2)
        {
          /* We only care about the first set of characters in each
             80-character row, so we don't need to parse the last few
             characters anyway. */
          fulllen=strlen(fullheader)-12;
          for(i=0;i<fulllen;++i)
            if( strncmp(fullheader+i, "CRVAL1  = '", 11) == 0 )
              fprintf(stderr, "WARNING: WCS Keyword values are not "
                      "numbers.\n\n"
                      "WARNING: The values to the WCS-related keywords are "
                      "enclosed in single-quotes. In the FITS standard "
                      "this is how string values are stored, therefore "
                      "WCSLIB is unable to read them AND WILL PUT ZERO IN "
                      "THEIR PLACE (creating a wrong WCS in the output). "
                      "Please update the respective keywords of the input "
                      "to be numbers (see next line).\n\n"
                      "WARNING: You can do this with Gnuastro's `astfits' "
                      "program and the `--update' option. The minimal WCS "
                      "keywords that need a numerical value are: `CRVAL1', "
                      "`CRVAL2', `CRPIX1', `CRPIX2', `EQUINOX' and "
                      "`CD%%_%%' (or `PC%%_%%', where the %% are integers), "
                      "please see the FITS standard, and inspect your FITS "
                      "file to identify the full set of keywords that you "
                      "need correct (for example PV%%_%% keywords).\n\n");
        }

      /* CTYPE is a mandatory WCS keyword, so if it hasn't been given (its
         '\0'), then the headers didn't have a WCS structure. However,
         WCSLIB still fills in the basic information (for example the
         dimensionality of the dataset). */
      if(outwcs->ctype[0][0]=='\0')
        {
          wcsfree(outwcs);
          outwcs=NULL;
          nwcs=0;
        }
      else
        {
          /* For a check.
          printf("flag: %d\n", outwcs->flag);
          printf("naxis: %d\n", outwcs->naxis);
          printf("crpix: %f, %f\n", outwcs->crpix[0], outwcs->crpix[1]);
          printf("pc: %f, %f, %f, %f\n", outwcs->pc[0], outwcs->pc[1], outwcs->pc[2],
                 outwcs->pc[3]);
          printf("cdelt: %f, %f\n", outwcs->cdelt[0], outwcs->cdelt[1]);
          printf("crval: %f, %f\n", outwcs->crval[0], outwcs->crval[1]);
          printf("cunit: %s, %s\n", outwcs->cunit[0], outwcs->cunit[1]);
          printf("ctype: %s, %s\n", outwcs->ctype[0], outwcs->ctype[1]);
          printf("lonpole: %f\n", outwcs->lonpole);
          printf("latpole: %f\n", outwcs->latpole);
          */

          /* Set the WCS structure. */
          status=wcsset(outwcs);
          if(status)
            {
              fprintf(stderr, "\n##################\n"
                      "WCSLIB Warning: wcsset ERROR %d: %s.\n"
                      "##################\n",
                      status, wcs_errmsg[status]);
              wcsfree(outwcs);
              outwcs=NULL;
              nwcs=0;
            }
          else
            /* A correctly useful WCS is present. When no PC matrix
               elements were present in the header, the default PC matrix
               (a unity matrix) is used. In this case WCSLIB doesn't set
               `altlin' (and gives it a value of 0). In Gnuastro, later on,
               we might need to know the type of the matrix used, so in
               such a case, we will set `altlin' to 1. */
            if(outwcs->altlin==0) outwcs->altlin=1;
        }
    }

  // /* For a check.
    wcsprt(outwcs);
  // */

  /* Clean up and return. */
  status=0;
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in freeing the memory used to "
                      "keep all the headers");
  return outwcs;

}





struct wcsprm *
gal_wcs_sip2tpv(struct wcsprm *inwcs,
                gal_data_t *in,
                char *infile,
                char *inhdu)
{
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */
  size_t i, fulllen;
  int nwcs, sumcheck;
  int nkeys=0, status=0;
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  struct wcsprm *outwcs=NULL;
  double cd[2][2]={0};
  double pv1[115]={0}, pv2[15]={0};

  char *fullheader=gal_wcs_add_pvkeywords(inwcs, in, pv1, pv2, cd, infile, inhdu, &nkeys);


  // printf("%s\n", fullheader);
  /* WCSlib function to parse the FITS headers. */
  status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, &nwcs, &outwcs);
  if( outwcs == NULL )
  {

    fprintf(stderr, "\n##################\n"
            "WCSLIB Warning: wcspih ERROR %d: %s.\n"
            "##################\n",
            status, wcs_errmsg[status]);
    outwcs=NULL; nwcs=0;
  }

  /* Set the internal structure: */
  if(outwcs)
    {
      /* It may happen that the WCS-related keyword values are stored as
         strings (they have single-quotes around them). In this case,
         WCSLIB will read the CRPIX and CRVAL values as zero. When this
         happens do a small check and abort, while informing the user about
         the problem. */
      sumcheck=0;
      for(i=0;i<outwcs->naxis;++i)
        {sumcheck += (outwcs->crval[i]==0.0f) + (outwcs->crpix[i]==0.0f);}
      if(sumcheck==outwcs->naxis*2)
        {
          /* We only care about the first set of characters in each
             80-character row, so we don't need to parse the last few
             characters anyway. */
          fulllen=strlen(fullheader)-12;
          for(i=0;i<fulllen;++i)
            if( strncmp(fullheader+i, "CRVAL1  = '", 11) == 0 )
              fprintf(stderr, "WARNING: WCS Keyword values are not "
                      "numbers.\n\n"
                      "WARNING: The values to the WCS-related keywords are "
                      "enclosed in single-quotes. In the FITS standard "
                      "this is how string values are stored, therefore "
                      "WCSLIB is unable to read them AND WILL PUT ZERO IN "
                      "THEIR PLACE (creating a wrong WCS in the output). "
                      "Please update the respective keywords of the input "
                      "to be numbers (see next line).\n\n"
                      "WARNING: You can do this with Gnuastro's `astfits' "
                      "program and the `--update' option. The minimal WCS "
                      "keywords that need a numerical value are: `CRVAL1', "
                      "`CRVAL2', `CRPIX1', `CRPIX2', `EQUINOX' and "
                      "`CD%%_%%' (or `PC%%_%%', where the %% are integers), "
                      "please see the FITS standard, and inspect your FITS "
                      "file to identify the full set of keywords that you "
                      "need correct (for example PV%%_%% keywords).\n\n");
        }

      /* CTYPE is a mandatory WCS keyword, so if it hasn't been given (its
         '\0'), then the headers didn't have a WCS structure. However,
         WCSLIB still fills in the basic information (for example the
         dimensionality of the dataset). */
      if(outwcs->ctype[0][0]=='\0')
        {
          wcsfree(outwcs);
          outwcs=NULL;
          nwcs=0;
        }
      else
        {
          /* For a check.
          printf("flag: %d\n", outwcs->flag);
          printf("naxis: %d\n", outwcs->naxis);
          printf("crpix: %f, %f\n", outwcs->crpix[0], outwcs->crpix[1]);
          printf("pc: %f, %f, %f, %f\n", outwcs->pc[0], outwcs->pc[1], outwcs->pc[2],
                 outwcs->pc[3]);
          printf("cdelt: %f, %f\n", outwcs->cdelt[0], outwcs->cdelt[1]);
          printf("crval: %f, %f\n", outwcs->crval[0], outwcs->crval[1]);
          printf("cunit: %s, %s\n", outwcs->cunit[0], outwcs->cunit[1]);
          printf("ctype: %s, %s\n", outwcs->ctype[0], outwcs->ctype[1]);
          printf("lonpole: %f\n", outwcs->lonpole);
          printf("latpole: %f\n", outwcs->latpole);
          */

          /* Set the WCS structure. */
          status=wcsset(outwcs);
          if(status)
            {
              fprintf(stderr, "\n##################\n"
                      "WCSLIB Warning: wcsset ERROR %d: %s.\n"
                      "##################\n",
                      status, wcs_errmsg[status]);
              wcsfree(outwcs);
              outwcs=NULL;
              nwcs=0;
            }
          else
            /* A correctly useful WCS is present. When no PC matrix
               elements were present in the header, the default PC matrix
               (a unity matrix) is used. In this case WCSLIB doesn't set
               `altlin' (and gives it a value of 0). In Gnuastro, later on,
               we might need to know the type of the matrix used, so in
               such a case, we will set `altlin' to 1. */
            if(outwcs->altlin==0) outwcs->altlin=1;
        }
    }

  /* For a check.
    wcsprt(outwcs);
  */

  /* Clean up and return. */
  status=0;
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in freeing the memory used to "
                      "keep all the headers");
  return outwcs;

}






int main(){


    gal_data_t *out=NULL;
    // struct linprm *lin=NULL;
    // struct wcsprm *tempwcs=NULL;
    double *dparm;
    struct wcsprm *wcs=NULL;
    struct wcsprm *outwcs=NULL;
    // struct disprm *dispre=NULL;
    // gal_fits_list_key_t *headers;
    char *outfile="./test-fits/test-sip.fits";
    char *infile="./test-fits/test-pv.fits", *inhdu="1";
    char *pvinfile="./test-fits/test-sip.fits";
    char *pvoutfile="./test-fits/sip2pv.fits";

    int nwcs;

    /* Temporary defined variables. Ignore after testing.
    double tpvx[8][8]={0};
    double tpvy[8][8]={0};
    struct disprm *dis=NULL;
    size_t a_order=0, b_order=0;
    double crpix[2]={0};
    size_t naxis[2]={0};
    double tpvu[8][8]={0};
    double tpvv[8][8]={0};
    double a_coeff[5][5]={0};
    double b_coeff[5][5]={0};
    double pv1[15] = {0};
    double pv2[15] = {0};
    double cd[2][2]={0};
    */

  // This block for pv2sip
  // /* Read wcs from fits file. */
  // wcs=gal_wcs_read(infile, inhdu, 0, 0, &nwcs);
  // gal_wcs_decompose_pc_cdelt(wcs);

  // /* Read the data of the input file. */
  // out=gal_fits_img_read(infile, inhdu, -1, 1);

  // outwcs=gal_wcs_tpv2sip(wcs, out, infile, inhdu);
  // wcsprt(outwcs);

  // /* Write the modified header into the fits file.  */
  // out->wcs=outwcs;
  // out->nwcs=1;

  // gal_fits_img_write(out, outfile, NULL, NULL);
  // gal_data_free(out);


  // wcs=gal_wcs_read(outfile, inhdu, 0, 0, &nwcs);
  // gal_wcs_decompose_pc_cdelt(wcs);







 // This block for sip2pv

  // gal_wcs_get_sipparams(wcs, cd, a_coeff, b_coeff);
  // gal_wcs_calc_sipeq( cd, pv1, pv2, outfile, inhdu);
  wcs=gal_wcs_read(pvinfile, inhdu, 0, 0, &nwcs);
  // gal_wcs_decompose_pc_cdelt(wcs);

  out=gal_fits_img_read(pvinfile, inhdu, -1, 1);

  outwcs=gal_wcs_sip2tpv(wcs, out, outfile, inhdu);

  /* Write the modified header into the fits file.  */
  out->wcs=outwcs;
  out->nwcs=1;

  gal_fits_img_write(out, pvoutfile, NULL, NULL);
  gal_data_free(out);



  return 0;
}
