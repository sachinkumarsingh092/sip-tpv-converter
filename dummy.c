#include <stdio.h>
#include <stdlib.h>
#include <gnuastro/wcs.h>

#include <wcslib/wcshdr.h>

#define IOBUFLEN 2880    /* size in bytes of each IO buffer (DONT CHANGE!) */

/* global variables */
 
#define FLEN_FILENAME 1025 /* max length of a filename  */
#define FLEN_KEYWORD   75  /* max length of a keyword (HIERARCH convention) */
#define FLEN_CARD      81  /* length of a FITS header card */
#define FLEN_VALUE     71  /* max length of a keyword value string */
#define FLEN_COMMENT   73  /* max length of a keyword comment string */
#define FLEN_ERRMSG    81  /* max length of a FITSIO error message */
#define FLEN_STATUS    31  /* max length of a FITSIO status text string */
 

char *buildheadersip(struct wcsprm *wcs, size_t num,
                     size_t i, size_t j, double value)
{
    size_t numkey=50, n;
    char *fullheader, fmt[50];
    sprintf(fmt, "%%-8s=%%71f");
    
    fullheader = malloc(num*80);

    /* Make keywords */
    for(n=0; n<num; n++)
    sprintf(fullheader+FLEN_CARD*n-1, fmt, "A_%d_%d", 4.445454);
    // sprintf(fullheader+FLEN_CARD-1, fmt, "KEYNAME2", 78.9894);
    printf("%s\n", fullheader);

    return fullheader;
}


struct wcsprm *gal_wcs_tpv2sip(struct wcsprm *inwcs){
    struct wcsprm *outwcs;
    char *fullheader=buildheadersip(inwcs, 0, 0 ,0, 99);
    int nkeys=0, status=0;
    int relax    = WCSHDR_all; /* Macro: use all informal WCS extensions. */
    int ctrl     = 0;          /* Don't report why a keyword wasn't used. */
    int nreject  = 0;          /* Number of keywords rejected for syntax. */
    int nwcs;

    /* WCSlib function to parse the FITS headers. */
    status=wcspih(fullheader, nkeys, relax, ctrl, &nreject, &nwcs, &outwcs);
    if(status)
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

  /* Clean up and return. */
  status=0;
  if (fits_free_memory(fullheader, &status) )
    gal_fits_io_error(status, "problem in freeing the memory used to "
                      "keep all the headers");
  return outwcs;

    return outwcs;

}

int main(){
    size_t numkey=50;
    char *fullheader, fmt[50];
    sprintf(fmt, "%%-8s=%%71f");
    
    fullheader = malloc(50*80);

    /* Make keywords */
    sprintf(fullheader, fmt, "KEYNAME", 4.445454);
    sprintf(fullheader+FLEN_CARD-1, fmt, "KEYNAME2", 78.9894);
    printf("%s\n", fullheader);

    return 0;
}