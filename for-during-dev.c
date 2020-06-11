#include <stdio.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>

#include <wcslib/wcslib.h>

/* High-level program to check library during development. */

int
main(int argc, char *argv[])
{
  int nwcs;
  gal_data_t *out;
  struct wcsprm *wcs;
  char *inhdu="1";
  char *infile="~/gsoc/gnuastro_dev/gnuastro-test-files/test-fits/test-sip.fits";
  

  wcs=gal_wcs_read(infile, inhdu, 0, 0, &nwcs);
  gal_wcs_decompose_pc_cdelt(wcs);  

  out=gal_fits_img_read(infile, inhdu, -1, 1);
  // wcsprt(wcs);
  /* The main function to check. */
  out->wcs=gal_wcs_distortion_convert(wcs,
				                             GAL_WCS_DISTORTION_TPV,
				                             out->dsize);

	out->nwcs=1;

  gal_fits_img_write(out, "testing-lib.fits", NULL, NULL);
  gal_data_free(out);
  /* Clean up and return! */
}
