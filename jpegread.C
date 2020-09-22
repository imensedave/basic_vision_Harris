#include <iostream>
#include <set>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <memory.h>
//#include <malloc.h> 
#include <ctype.h>
#include <string.h>
#include <string>
#include <strings.h>
#include <sys/types.h> 



extern "C" {
   #include "jpeglib.h"
   #include <setjmp.h>
}

#include "Imgx.h"


// basic jpeg reader from libjpeg.

extern int jpeg_read( const char * filename, CxImage &im)
{

  
  FILE *fp = fopen(filename, "rb");
  if (fp==0) return( 0 );
  
  struct jpeg_decompress_struct jds;
  struct jpeg_error_mgr jerr;
  jds.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&jds);
  jpeg_stdio_src(&jds, fp);


  jpeg_read_header(&jds, true);
  if (jds.out_color_space != JCS_GRAYSCALE) jds.out_color_space = JCS_RGB;
  jds.output_components  = 3;
  jpeg_calc_output_dimensions(&jds);

  CxImage image(jds.output_height, jds.output_width);

  jpeg_start_decompress(&jds);
  JSAMPROW scan = new JSAMPLE [ 3 * jds.output_width ];
  for(unsigned int i = 0; i < jds.output_height; i++) {
    jpeg_read_scanlines(&jds, &scan, 1);
    JSAMPLE * scan_ptr = scan;

    short *img_ptrR = image.red()[i];
    short *img_ptrG = image.green()[i];
    short *img_ptrB = image.blue()[i];

    if (jds.out_color_space == JCS_GRAYSCALE) {
      for(unsigned int j = 0; j < jds.output_width; j++) {
	*img_ptrR++ = *scan_ptr;
	*img_ptrG++ = *scan_ptr;
	*img_ptrB++ = *scan_ptr;
	scan_ptr++;
      }
    }
    else{
      for(unsigned int j = 0; j < jds.output_width; j++) {
	*img_ptrR++ = *scan_ptr;  scan_ptr++;
	*img_ptrG++ = *scan_ptr;  scan_ptr++;
	*img_ptrB++ = *scan_ptr;  scan_ptr++;
      }
    }
  }

  jpeg_finish_decompress(&jds);
  jpeg_destroy_decompress(&jds);

  im = image;
  delete [] scan;
  fclose(fp);
  return( 1 );

}

extern int jpeg_read_grey( const char * filename, Matrix<short> &im)
{

  
  FILE *fp = fopen(filename, "rb");
  if (fp==0) return( 0 );
  
  struct jpeg_decompress_struct jds;
  struct jpeg_error_mgr jerr;
  jds.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&jds);
  jpeg_stdio_src(&jds, fp);


  jpeg_read_header(&jds, true);
  if (jds.out_color_space != JCS_GRAYSCALE) jds.out_color_space = JCS_RGB;
  jds.output_components  = 3;
  jpeg_calc_output_dimensions(&jds);

  //CxImage image(jds.output_height, jds.output_width);
  im.reshape(jds.output_height, jds.output_width);

  jpeg_start_decompress(&jds);
  JSAMPROW scan = new JSAMPLE [ 3 * jds.output_width ];
  for(unsigned int i = 0; i < jds.output_height; i++) {
    jpeg_read_scanlines(&jds, &scan, 1);
    JSAMPLE * scan_ptr = scan;

    short *pim = im[i];

    if (jds.out_color_space == JCS_GRAYSCALE) {
      for(unsigned int j = 0; j < jds.output_width; j++) {
	*pim++ = *scan_ptr;
	scan_ptr++;
      }
    }
    else{
      for(unsigned int j = 0; j < jds.output_width; j++) {
	//tmp = (float) 0.4* ( max( *R, *B)) +  0.6* (*G);// +  0.114* (*B) + 0.5; // matlabs version.
	int r =  *scan_ptr; scan_ptr++;
	int g =  *scan_ptr; scan_ptr++;
	int b =  *scan_ptr; scan_ptr++;

	
	int val = (52*(max(r,b)) + 76*g)>>7;
	*pim = (short) val;
	pim++;
      }
    }
  }

  jpeg_finish_decompress(&jds);
  jpeg_destroy_decompress(&jds);

  delete [] scan;
  fclose(fp);
  return( 1 );

}

/* from the jpeg library example routine
 *
 *
 * 60 looks like good quality
 *
 */


extern int jpeg_write( const char * filename, int quality, CxImage &cim)
{
   
  /* This struct contains the JPEG compression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   * It is possible to have several such structures, representing multiple
   * compression/decompression processes, in existence at once.  We refer
   * to any one struct (and its associated working data) as a "JPEG object".
   */
  struct jpeg_compress_struct cinfo;
  /* This struct represents a JPEG error handler.  It is declared separately
   * because applications often want to supply a specialized error handler
   * (see the second half of this file for an example).  But here we just
   * take the easy way out and use the standard error handler, which will
   * print a message on stderr and call exit() if compression fails.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct jpeg_error_mgr jerr;
  /* More stuff */
  FILE * outfile;		/* target file */
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */

  /* Step 1: allocate and initialize JPEG compression object */

  /* We have to set up the error handler first, in case the initialization
   * step fails.  (Unlikely, but it could happen if you are out of memory.)
   * This routine fills in the contents of struct jerr, and returns jerr's
   * address which we place into the link field in cinfo.
   */
  cinfo.err = jpeg_std_error(&jerr);
  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&cinfo);

  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */

  /* Here we use the library-supplied code to send compressed data to a
   * stdio stream.  You can also write your own code to do something else.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to write binary files.
   */
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }
  jpeg_stdio_dest(&cinfo, outfile);

  /* Step 3: set parameters for compression */

  /* First we supply a description of the input image.
   * Four fields of the cinfo struct must be filled in:
   */
  cinfo.image_width = cim.cols(); //image_width; 	/* image width and height, in pixels */
  cinfo.image_height = cim.rows(); //image_height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&cinfo);
  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Step 4: Start compressor */

  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&cinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
  row_stride = cim.cols() * 3;	/* JSAMPLEs per row in image_buffer */
  short *pr,*pg,*pb;
  int r=0, nc = cim.cols();

  JSAMPROW scan = new JSAMPLE [ 3 * cim.cols() ];

  while (cinfo.next_scanline < cinfo.image_height) {
    JSAMPLE * scan_ptr = scan;
    pr = cim.red()[r];
    pg = cim.green()[r];
    pb = cim.blue()[r];
    r++;
    for( int q=0; q<nc; q++, pr++, pg++, pb++){
      *scan_ptr = *pr; scan_ptr++;
      *scan_ptr = *pg; scan_ptr++;
      *scan_ptr = *pb; scan_ptr++;
    }

    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    //    row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
    scan_ptr = scan;
    (void) jpeg_write_scanlines(&cinfo, &scan_ptr, 1);
  }

  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);
  /* After finish_compress, we can close the output file. */
  fclose(outfile);

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_compress(&cinfo);

  /* And we're done! */
}
