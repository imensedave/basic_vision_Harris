#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <memory.h>
#include <ctype.h>
#include <string>
#include <sys/types.h>
#include <sys/time.h>

#include <jpeglib.h>

#include <ctime>

#include "Imgx.h"
#include "lodepng.h"
#include "bmp.h"
#include "IM_struct.h"
#include "IM_func.h"


/* copyright 2020 Imense Ltds
 * 
 * opensource image/ vision example project.
 * the code in this file is released under the apache2 license.
 

*/

using namespace std;


int main(int argc, char *argv[])
{

  if( argc == 1 ){
    cerr << "USAGE imx <file.jpg>" << endl;
    cerr << "read in various image file formats" << endl;
    cerr << "perform something a bit like Canny edge detection " << endl;
    cerr << "write something to a file you can view in mutlab " << endl;

    exit(2);
  }


  CxImage cim(1,1); // colour image
  int nr, nc , area, goat_balls = 1;
  unsigned char *duff;
  Matrix<short> tim;
  Matrix<short> grey;


  // sometimes you want to time things
  struct timeval t0, t1,t2; gettimeofday(&t0,NULL);

  // use input filename ending to determine file type.
  std::string fname(argv[1]), end1(".jpg"), end2(".JPG"), end3(".bmp"), 
    end4(".BMP"), end7("jpeg"), end8(".png"), end12(".pgm");

  
  if( fname.find(end1) != fname.npos || fname.find(end2) != fname.npos || fname.find(end7) != fname.npos ){
    
    // read to a colour image:
    int foo  = jpeg_read( argv[1], cim);
    if( foo == 0 ){
      cerr << "jpeg read error " << argv[1] << endl;
      exit(2);
    }
    
    // read directly to a greyscale image
    int foo2  = jpeg_read_grey( argv[1], grey);
    if( foo2 == 0 ){
      cerr << "jpeg read error " << argv[1] << endl;
      exit(2);
    }
    nr = grey.rows(); nc = grey.cols(); area = nr*nc;
  }

  else if( fname.find(end3) != fname.npos || fname.find(end4) != fname.npos ){
    BITMAPINFOHEADER2 bitmapInfoHeader;
    duff = LoadBitmapFile(argv[1], &bitmapInfoHeader);
    if( duff == NULL ){
      cerr << "bmp read error " << argv[1] << endl;
      exit(2);
    }
    nr = bitmapInfoHeader.biHeight;
    nc = bitmapInfoHeader.biWidth;
    goat_balls = 0;
    grey.reshape(nr,nc);
  }
  else if( fname.find(end12) != fname.npos ){
    int ham = read_pgm( tim, argv[1]);
    if( ham <= 0 ){
      cerr << "pgm read error " << argv[1] << endl;
      exit(2);
    }
    nr = tim.rows();
    nc = tim.cols();
    goat_balls = 3;
    grey.reshape(nr,nc);
  }
  else if(fname.find(end8) != fname.npos ){
    int knee = decodeTwoSteps(argv[1], cim);
    if( knee <= 0 ){
      cerr << "png read error" << argv[1] << endl;
      exit(2);
    }
    nr = cim.rows();
    nc = cim.cols();
    goat_balls = 2;
    grey.reshape(nr,nc);
  }
  else{
    cerr << "filetype not supported " << argv[1]  << endl;
    exit(2);
  }

  if( goat_balls == 2 ){
    // create a greyscale image from the colout image.
    C_to_BW( cim, grey);
  }
  else if( goat_balls == 3) {
    grey = tim;
  }
  else if(goat_balls == 0 ){
    unsigned char  *p=duff;
    grey.reshape(nr,nc);
    for( int r=nr-1; r>=0; r--){
      short *p2 = grey[r], *pend = p2+nc;
      for( ; p2<pend; p++, p2++){
	(*p2) = (int) (*p);
      }
    }
  }

  //  grey.write_rlm_ascii("grey.rlm");
  Matrix<short> edges(nr,nc);
  
  // routine to perform canny type edge detection.
  std::vector<PT> pts;
  float sigma = 2.5, scale = 0.04; 
  Harris_corners(grey, pts, sigma, scale );
  
    
  gettimeofday(&t1,NULL);
  float diff = t1.tv_usec + 1000000*t1.tv_sec -( t0.tv_usec + 1000000*t0.tv_sec);

  cerr << "processing time in millie seconds: " << diff/1000 << endl << endl;
  
  write_corners(pts, "corners.mat");

  
  exit(0);
}

