
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <set>
#include <string>
#include <cfloat>
#include "Imgx.h"
#include "IM_struct.h"
#include "IM_func.h"



// copyright David Sinclair 2020


/* released under the Apache2 license.

   https://www.apache.org/licenses/LICENSE-2.0


   dr.sinclair@gmail.com
*/


void  make_gauss(float *mask, int num, float sigma);
void  conv1r( Matrix<short> &im, Matrix<short> &out, float *mask, int num );
void  conv1c( Matrix<short> &im, Matrix<short> &out, float *mask, int num );




#define debug_flag 0

/* Harris and Stevens style corner detector with scale given by size of Gaussian convolution kernel.

   http://www.bmva.org/bmvc/1988/avc-88-023.pdf
 
   written for ease of understanding of code not blinding speed.


*/

void Harris_corners( Matrix<short> &im, std::vector<PT> &pts, float sigma, float scal)
{
  int nr = im.rows(), nc = im.cols();
  int nr2 = nr-3, nc2 = nc-3;
  Matrix<short> sim(im), tmp(im);

  float mask[7];
  make_gauss(mask, 7, sigma);

  Matrix<short> drx(nr,nc), dcx(nr,nc), drcx(nr,nc);
  
  for (int c=3; c<nc2; c++){
    for( int r=3; r<nr2; r++){

      int dr = 2*(sim[r+1][c] - sim[r-1][c])
	+ sim[r+1][c+1] - sim[r-1][c-1]
	+ sim[r+1][c+1] - sim[r-1][c+1];
      int dc = 2*(sim[r][c+1] - sim[r][c-1]) + 
	sim[r-1][c+1] - sim[r-1][c-1] +
	sim[r+1][c+1] - sim[r+1][c-1];

      dr >>= 2;
      dc >>= 2;
      
      drx[r][c] = dr*dr;
      dcx[r][c] = dc*dc;
      drcx[r][c] = dr*dc;
    }
  }

  // smooth image by convolution with Gaussian mask.
  
  conv1r( drx, tmp, mask, 7 );
  conv1c( tmp, drx, mask, 7 );

  conv1r( dcx, tmp, mask, 7 );
  conv1c( tmp, dcx, mask, 7 );
  
  conv1r( drcx, tmp, mask, 7 );
  conv1c( tmp, drcx, mask, 7 );

  Matrix<int> Hmat(nr,nc);

  nc2 - nc-8; nr2 = nr-8;
  for (int c=8; c<nc2; c++){
    for( int r=8; r<nr2; r++){

      int A = drx[r][c];
      int B = dcx[r][c];
      int C = drcx[r][c];
      
      float H = A*B -C*C - scal* (A+B)*(A+B);

      Hmat[r][c] = (short) H;

    }
  }

  if( debug_flag ){
    drx.write_rlm_ascii("dr.rlm");
    dcx.write_rlm_ascii("dc.rlm");
    drcx.write_rlm_ascii("drc.rlm");
    tmp.write_rlm_ascii("harris.rlm");
  }
  
  std::vector<PT> corns;
  float v;

  // Yeah Olde n-max suppression (8 connected).

  for(int r=2; r<nr-2; r++){
    for(int c=2; c<nc-2; c++){
      v =Hmat[r][c];
      // if( v > cmax ){
      if( v > 10000 ){
	if( Hmat[ r-1][c+1] <= v)
	  if( Hmat[ r-1][c  ] <= v)
            if( Hmat[ r-1][c-1] <= v)
	      if( Hmat[ r  ][c-1] <= v)
		
		if( Hmat[ r+1][c-1] < v)
		  if( Hmat[ r+1][c  ] < v)
		    if( Hmat[ r+1][c+1] < v)
		      if( Hmat[ r  ][c+1] < v){
			PT t(r,c, v);
			corns.push_back(t);
		      }
      }
    }
  }
    
  sort( corns.begin(), corns.end());
  int numc = corns.size();
  if( numc > 2000){
    numc = 2000;
  }
  
  for( int q=0; q<numc; q++){
    corns[q].v = sqrt(sqrt(corns[q].v));
    pts.push_back(corns[q]);
  }

  if( debug_flag > 0 ){
    cerr << numc << " Harris corners found " << endl;
  }

  






  
  return;
}







/* convolve up and down a coloumn along rows
   (think snails reading a book),
 
   mask is a half gaussian.  ^\_

*/

void  conv1r(Matrix<short> &im, Matrix<short> &out, float *mask, int num )
{
   int nr = im.rows(), nc = im.cols();
   int nc2 = nc-num;

   for( int r=0; r<nr; r++){
     for (int c=num; c<nc2; c++){

       float v = (float) mask[0]* im[r][c];
       for( int k=1; k<num; k++){
	 v += (im[r][c+k] + im[r][c-k])* mask[k];
       }
       out[r][c] = short( v+0.0499);
     }
   }
   return;
}

/* convolve up and down a coloumn.
   like gravey dribbling down a pinstriped shirt.
 
   mask is a half gaussian.  ^\_

*/

void  conv1c(Matrix<short> &im, Matrix<short> &out, float *mask, int num )
{
   int nr = im.rows(), nc = im.cols();
   int nr2 = nr-num;

   for (int c=0; c<nc; c++){
     for( int r=num; r<nr2; r++){

       float v = mask[0]* im[r][c];
       for( int k=1; k<num; k++){
	 v += (im[r+k][c] + im[r-k][c])* mask[k];
       }
       out[r][c] = short( v+0.0499);
     }
   }
   return;
}


/* mask is a half Gaussian with value at 0 at mask[0].
 * 

 */

void make_gauss(float *mask, int num, float sigma)
{
   int nx = num/2;
   
   memset((void *) mask, 0, sizeof(float)*num);
   
   float total;
   float sigma2 = 2*sigma*sigma;

   float *mx = mask;
   
   for(int i=0; i<=num; i++ ){
      float tmpx = i - 0.5;
      for(int j=0; j<11; j++){
         mx[i] += exp( - tmpx * tmpx/ (sigma2) ) ;
         tmpx += 0.1;
      }
      mx[i] = mx[i]/10.0;
   }
   total = mx[0];
   for(int i=1; i<=nx; i++ )
      total += 2.0*mx[i];

   for(int i=0; i<num; i++ ) {
      mx[i] /= total;
   }
   
   return;
}

void write_corners(std::vector<PT> &ts, char * fname){
	std::ofstream out(fname, ios::out);
	
	int nt = ts.size();
	out << nt << " 3" << std::endl;
	
	for (int i = 0; i < nt; i++){
		out << ts[i].r << " " <<  ts[i].c << " " << ts[i].v << std::endl;
	}
	
	out.close();
	
	return;
};
