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

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/permutation.h>





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
gal_wcs_get_params(struct wcsprm *wcs, double cd[2][2], \
                   double *pv1, double *pv2)
{
  size_t i, j,k,index=0;
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
  tpvy[0][0] = pv2[0];
  tpvy[0][1] = pv2[1];
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
gal_wcs_real_tpveq(double cd[2][2], double tpvu[8][8], double tpvv[8][8],
                   char *infile, char *inhdu)
{
  /* tpvu and tpvv are u-v distortion equations in TPV convention. */
  int a;
  size_t i,j;
  double determinant=0;
  struct wcsprm *wcs=NULL;
  double cd_inv[2][2]={0};
  double pv1[GAL_WCS_MAX_PVSIZE] = {0};
  double pv2[GAL_WCS_MAX_PVSIZE] = {0};
  double tpv1[8][8]={0}, tpv2[8][8]={0};


  wcs=gal_wcs_read(infile, inhdu, 0, 0, &a); 

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

                  |tpvu| = cd_inv*|tpv1|
                  |tpvv|          |tpv2|
    though intermidate polynomial equations have to be calculated prior
    to multiplycation with cd_inv.
      */

  for(i=0; i<8; i++)
    for(j=0; j<8; j++){
      tpvu[i][j]=cd_inv[0][0]*tpv1[i][j]+cd_inv[0][1]*tpv2[i][j];
      tpvv[i][j]=cd_inv[1][0]*tpv1[i][j]+cd_inv[1][1]*tpv2[i][j];

      /*For a check:
      printf("%.10f %.10f\n", tpv1[i][j], tpv2[i][j]);
      printf("%.10lf, %.10lf\n", tpvu[i][j], tpvv[i][j]);
      */
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






char *
gal_wcs_add_revkeywords(){
  /*TODO*/
}







char *  
gal_wcs_add_sipkeywords(struct wcsprm *wcs, double tpvu[8][8], double tpvv[8][8], 
                        char *infile, char *inhdu, int add_reverse, int *nkeys)
{
  double val=0;
  uint8_t i, j, k=0;
  double cd[2][2] = {0};
  int size = wcs->naxis;
  size_t a_order=0, b_order=0;
  size_t m, n, num=0, numkey=100;
  char *fullheader, fmt[50], sipkey[8], keyaxis[9], pcaxis[10];
  
  *nkeys = 0;

  /* The format for each card. */
  sprintf(fmt, "%%-8s= %%20.12E%%50s");

  gal_wcs_real_tpveq(cd, tpvu, tpvv, infile, inhdu);
  
  
  /* Allcate memory for cards. */
  fullheader = malloc(numkey*80);

  /* Add other necessary cards. */
  sprintf(fullheader+(FLEN_CARD-1)*num++, "%-8s= %20d%50s", "WCSAXES", wcs->naxis, "");
  

  /* Include a loop for each axes. */
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


  for(m=0; m<8; m++)
    for(n=0; n<8; n++)
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
    gal_wcs_add_revkeywords();


  *nkeys = num;

  /*For a check.
  printf("%s\n", fullheader);
  */

  return fullheader;

}





struct wcsprm *
gal_wcs_tpv2sip(struct wcsprm *inwcs, char *infile, char *inhdu)
{
  int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
  int nreject  = 0;          /* Number of keywords rejected for syntax. */
  size_t i, fulllen;
  int nwcs, sumcheck;
  int nkeys=0, status=0;
  int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
  struct wcsprm *outwcs=NULL;
  double tpvu[8][8] ={0}, tpvv[8][8]={0};
  
  char *fullheader=gal_wcs_add_sipkeywords(inwcs, tpvu, tpvv, infile, inhdu, 0, &nkeys);


  printf("%s\n", fullheader);
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

  /*For a check.
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
    char *wcsstr;
    gal_data_t *out=NULL;
    // struct linprm *lin=NULL;
    struct wcsprm *wcs=NULL;
    struct wcsprm *outwcs=NULL;
    // struct disprm *disseq=NULL, *tdisp=NULL;
    // gal_fits_list_key_t *headers;
    char *outfile="test-sip.fits";
    char *infile="test-pv.fits", *inhdu="1";

    int nwcs, nkeyrec, dstatus, wstatus;

    /* Temporary defined variables. Ignored while testing.
    double cd[2][2]={0};
    double tpvx[8][8]={0};
    double tpvy[8][8]={0};
    double tpvu[8][8]={0};
    double tpvv[8][8]={0};
    struct disprm *dis=NULL;
    double pv1[GAL_WCS_MAX_PVSIZE] = {0};
    double pv2[GAL_WCS_MAX_PVSIZE] = {0};
    */

  /* Read wcs from fits file. */
  wcs=gal_wcs_read(infile, inhdu, 0, 0, &nwcs);
  gal_wcs_decompose_pc_cdelt(wcs);


  // gal_wcs_real_tpveq(cd, tpvu, tpvv, infile, inhdu);
  // gal_wcs_add_sipkeywords(tpvu, tpvv, 0);

  outwcs=gal_wcs_tpv2sip(wcs, infile, inhdu);
  wcsprt(outwcs);

  /* Read the data of the input file. */
  out=gal_fits_img_read(infile, inhdu, -1, 1);

  /* Write the modified header into the fits file.  */
  out->wcs=outwcs;
  out->nwcs=1;

  gal_fits_img_write(out, outfile, NULL, NULL);
  // gal_data_free(out);

  
  return 0;
}
