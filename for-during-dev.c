#include <stdio.h>
#include <stdlib.h>

#include <gnaustro/wcs.h>
#include <gnaustro/fits.h>

/* High-level program to check library during development. */

int
main(int argc, char *argv[])
{
  gal_data_t *out;
  struct wcsprm *wcs;
  char *inhdu=argv[2];
  char *infile=argv[1];

  wcs=gal_wcs_read(infile, inhdu, 0, 0, &nwcs);

  out=gal_fits_img_read(infile, inhdu, -1, 1);

  /* The main function to check. */
  out->wcs=gal_wcs_distortion_convert(wcs,
				      GAL_WCS_DISTORTION_SIP,
				      out->dsize);

  gal_fits_img_write(out, pvoutfile, NULL, NULL);

  /* Clean up and return! */
}
