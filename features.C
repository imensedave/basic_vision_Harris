#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector.h>
#include <set.h>
#include <string.h>
#include <string>
#include <stack.h>
#include <algo.h>
#include <time.h>


#include "structures.h"
#include "Vimage.h"
#include "utils.h"
#include "Lines.h"


/* corner related routines.
 *
 * copyright 2008 Imense Ltd 
 */


/* 
 * plessy style corners 
 *
 * somewhat memory intensive (could be handled as a local operator
 * but isn't).
 *
 */
void  cornerx(Matrix<float> &im, std::vector<PT> &corns, float sigma, float scal)
{
   int nr = im.rows(), nc = im.cols();
   Matrix<float> dx(nr,nc), dy(nr,nc), dxy(nr,nc), tim(nr,nc);
   float grad[3]={-1,0,1};
   float *p, *p1, *p2, *p3;
   
   conv1r(im, dx, grad, 3 );
   //write_mat( dx, "dx.mut" );
   
   conv1c(im, dy, grad, 3 );
   //write_mat( dy, "dy.mut" );
   
   for( int i=0; i<nr; i++){
      dx[i][0] = 0;      dy[i][0] = 0;
      dx[i][nc-1] = 0;      dy[i][nc-1] = 0;
   }
   for( int i=0; i<nc; i++){
      dx[0][i] = 0;      dy[0][i] = 0;
      dx[nr-1][i] = 0;      dy[nr-1][i] = 0;
   }
      
   
   
   // make dxy
   int nrc = nr*nc;
   p = dxy[0]; p1 = dy[0]; p2 = dx[0];
   float *pend=p+nrc;
   for(; p<pend; p++, p1++, p2++){
      *p = (*p1)* (*p2);
   }
   //write_mat( dxy, "dxy.mut" );
   
   float mask[7];
   make_gauss(mask, 7, sigma);

   conv1r(dxy,tim,mask,7);
   conv1c(tim,dxy,mask,7);
   
   
   dx.square();
   conv1r(dx,tim,mask,7);
   conv1c(tim,dx,mask,7);
   dy.square();
   conv1r(dy,tim,mask,7);
   conv1c(tim,dy,mask,7);

   // compute cornerness! 
   p = tim[0];
   p3 = dxy[0]; p1 = dx[0]; p2 = dy[0];
   pend=p+nrc;
   for(; p<pend; p++, p1++, p2++, p3++){
      *p = (*p1)* (*p2) - (*p3)* (*p3) - (*p1+*p2) * scal ;
   }
   //write_mat( tim, "corns.mut" );
   
   for( int i=0; i<nr; i++){
      tim[i][0] = 0;      tim[i][1] = 0;
      tim[i][nc-1] = 0;   tim[i][nc-2] = 0;
   }
   for( int i=0; i<nc; i++){
      tim[0][i] = 0;      tim[1][i] = 0;
      tim[nr-1][i] = 0;   tim[nr-2][i] = 0;
   }
   
   //write_mat( tim, "corns.mut" );
   
   /*   p = tim[0]; 
   pend=p+nrc;
   float cmax = -1;
   float cmin = 0.0001;
   for(; p<pend; p++){
      if( *p > cmax)
         cmax = *p;
   }
   
   
   cmax *=cmin; */
   float v;
   for(int r=2; r<nr-2; r++){
      for(int c=2; c<nc-2; c++){
         v =tim[r][c];
         // if( v > cmax ){
         if( v > 10000 ){
            if( tim[ r-1][c+1] <= v)
            if( tim[ r-1][c  ] <= v)
            if( tim[ r-1][c-1] <= v)
            if( tim[ r  ][c-1] <= v)

            if( tim[ r+1][c-1] < v)
            if( tim[ r+1][c  ] < v)
            if( tim[ r+1][c+1] < v)
            if( tim[ r  ][c+1] < v){
	      PT t(r,c, v);
               corns.push_back(t);
            }
         }
      }
   }
    
   sort( corns.begin(), corns.end());
   int numc = corns.size();
   if( numc > 2000){
     std::vector<PT> cornx;
     for( int q=0; q<2000; q++){
       corns[q].v = sqrt(sqrt(corns[q].v));
       cornx.push_back(corns[q]);
     }
     corns = cornx;
     fprintf(stderr, "n0 corners %d,  val_1000 = %.2f\n", numc, corns[1999].v);
   }
   else
     fprintf(stderr, "n0 corners %d,  val_xxx = %.2f\n", numc, corns[numc-1].v);

   

        
   return;
}

/* 
 * out is assumed to be of the right size on entry
 * 
 * mask must be of odd size
 */

// convolve along row.
void  conv1r(Matrix<float> &im, Matrix<float> &out, float *mask, int num )
{
   int nr = im.rows(), nc = im.cols();
   int ox = num/2;
   int nc2 = nc+num*2;
 
   float tmp[nc2], *f;
   float *mx = mask+ox;
   
	for (int i=0; i<nr;i++){
      float *p = im[i], *pend = p+nc;

      memset((void *) tmp, 0, sizeof(float)*nc2);
      f = tmp+num;
      for(p = im[i]; p<pend; p++, f++){
         for(int x=-ox; x<=ox; x++){
            f[-x] += mx[x] *(*p);
			}
      }

      /* fix the boundary */
      f = tmp+num;
      for(int x=0; x<=ox; x++){
         f[x] += f[-x-1];
      }
      for(int x=0; x<=ox; x++){
         f[nc-1-x] += f[nc+x];
      }
      // copy buffer.
      p = out[i]; 
      f = tmp+num;
      pend = p+nc;
      for(; p<pend; p++, f++){
         *p = *f;
      }
	}
		
	return;
}


// convolve up and down a coloumn.
void  conv1c(Matrix<float> &im, Matrix<float> &out, float *mask, int num )
{
   int nr = im.rows(), nc = im.cols();
   int ox = num/2, nrc = nr*nc;
   int nr2 = nr+num*2;
   float *p, *pend;
    
   float tmp[nr2], *f;
   float *mx = mask+ox;
   
	for (int i=0; i<nc;i++){

      memset((void *) tmp, 0, sizeof(float)*nr2);
      f = tmp+num;
      p = im[0]+i; 
      pend = p+nrc;
      for(; p<pend; p+=nc, f++){
         for(int x=-ox; x<=ox; x++){
            f[-x] += mx[x] *(*p);
			}
      }

      /* fix the boundary */
      f = tmp+num;
      for(int x=0; x<=ox; x++){
         f[x] += f[-x-1];
      }
      for(int x=0; x<=ox; x++){
         f[nr-1-x] += f[nr+x];
      }
      // copy buffer.
      p = out[0]+i; 
      f = tmp+num;
      pend = p+nrc;
      for(; p<pend; p+=nc, f++){
         *p = *f;
      }
	}
		
	return;
}


/*
 * num is assumed to be an odd number
 *
 *
 */

void make_gauss(float *mask, int num, float sigma)
{
   int nx = num/2;
   float *mx =mask+nx;
   
   memset((void *) mask, 0, sizeof(float)*num);
   
   float total;
   float sigma2 = 2*sigma*sigma;

   for(int i=0; i<=nx; i++ ){
      float tmpx = i - 0.5;
      for(int j=0; j<11; j++){
         mx[i] += exp( - tmpx * tmpx/ (sigma2) ) ;
         tmpx += 0.1;
      }
      mx[i] = mx[i]/10.0;
      //fprintf(stderr, " %0.6f\n", mx[i]);
   }
   total = mx[0];
   for(int i=1; i<=nx; i++ )
      total += 2.0*mx[i];

   mx[0] /= total;
   //cerr << mx[0] << " " ;
   for(int i=1; i<=nx; i++ ) {
      mx[i] /= total;
      mx[-i] = mx[i];
     // cerr << mx[i] << " " ;
   }
//cerr << endl;
   
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

void write_mat( Matrix <float> &m, char * fname )
{
   int nr = m.rows(), nc = m.cols();
   FILE *fp = fopen(fname, "w");
    
   fprintf(fp,"%d %d\n",nr, nc);

   for (int i=0; i<nr; i++){
      for (int j=0; j<nc; j++){
	if( m[i][j] > 0.01)
	  fprintf(fp," %.2f", m[i][j]);
	else
	  fprintf(fp," %d", 0);   // tuned to harris corner responce.
      }
      fprintf(fp,"\n");
   }

   fclose(fp);
   
   return;
}

void write_mat( Matrix <int> &m, char * fname )
{
  int nr = m.rows(), nc = m.cols();
  FILE *fp = fopen(fname, "w");
  
  fprintf(fp,"%d %d\n",nr, nc);
  
  for (int i=0; i<nr; i++){
    for (int j=0; j<nc; j++){
      fprintf(fp," %d", m[i][j]);
    }
    fprintf(fp,"\n");
  }
  
  fclose(fp);
  
  return;
}

void C_to_float( CImage &cim, Matrix<float> &out )
{

  int nr=cim.rows();
  int nc=cim.cols();
  if( ( nr != (int) out.rows()) || ( nc != (int) out.cols() ) ){
    Matrix<float> df(nr, nc );
    out = df;
  }
  
  int *R, *G, *B;
  float *ptr, *pend;
  ptr = out[0]; pend = ptr + nr*nc;
  R = cim.red()[0];
  G = cim.green()[0];
  B = cim.blue()[0];
  

  for (; ptr<pend; ptr++, R++, G++, B++){
	  //tmp = (float) 0.34375* (*R) +  0.5* (*G) +  0.15625* (*B); 
	  *ptr = (float) 0.2989* (*R) +  0.587* (*G) +  0.114* (*B) ; // matlabs version.
  }
}



/* start from grey image to keep things simple
 * one scale edge detection works well enough.
 *
 *
 * returns chains and edge map and turning angle map.
 * 
 * edges are given to sub pixel accuracy (quadratic fit) to reduce alyassing woes.
 *
 * image is assumed to be lightly smoothen before start.

 */

int canny_tangle(Matrix<float> &im, float low, float high, BWImage &edge, Matrix<float> &tangle,
		 std::vector<Exf_chain> &chains, Matrix<float> &drx, Matrix<float> &dcx, 
		 Matrix<float> &magx, std::vector<Exf> &hc_pts, std::vector<Line> &lines )
{
  int i,j, *ptr, *ptc;
  float mag, *pend, *r, *ptf, dr,dc;
  chains.clear();
		
  int nr = im.rows(); int br = nr-1;
  int nc = im.cols(); int bc = nc-1;
		
  /*  Matrix<float> magx(nr,nc);
  Matrix<float> drx(nr,nc);
  Matrix<float> dcx(nr,nc); */
  Matrix<float> smot(nr,nc);
  

  // compute derivatives
  for(i=1; i<br; i++){
    r = im[i];
    pend = r + bc;
    r += 1;
    
    for( j=1; r<pend; r++, j++ ){
      dcx[i][j] = (*(r+1 )) - (*(r-1 ));
      drx[i][j] = (*(r+nc)) - (*(r-nc));
      magx[i][j] = dcx[i][j]*dcx[i][j] + drx[i][j]*drx[i][j];
    }
  }
  //write_f_mat2(magx, "edge_mag");
  //write_f_mat2(drx, "drx");
  //write_f_mat2(dcx, "dcx");
		
  edge.clear();
  smot.clear();
	
  // non-maximum supression
  int count;
  float val, theta;
  //  float p1 = -7.0*3.141592654/8.0, p2=-5.0*3.141592654/8.0, p3=-3*3.141592654/8.0, p4=-3.141592654/8.0, p5=3.141592654/8.0, p6=3.0*3.141592654/8.0, p7=5.0*3.141592654/8.0, p8=7.0*3.141592654/8.0;
  float p1 = -6.8*3.141592654/8.0, p2=-5.2*3.141592654/8.0, p3=-2.8*3.141592654/8.0, p4=-1.2*3.141592654/8.0, 
        p5 =  1.2*3.141592654/8.0, p6=2.8*3.141592654/8.0,  p7= 5.2*3.141592654/8.0, p8= 6.8*3.141592654/8.0;
  float nms_strength=0, vr, vc, x, r2 = 1.4142;
  vector<int> vi, vj;
	
  count = 0;
  for(i=1; i<br; i++){
    ptf = magx[i];
    pend  = ptf + bc;
    ptf  += 1;
    
    for(int j=1 ; ptf<pend; ptf++, j++){
      if( *ptf > low ){
	vr = (float) drx[i][j];
	vc = (float) dcx[i][j];
	val = *ptf;
	
	theta = atan2(vr,vc);
	if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	  if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	  if( magx[i-1][j-1]<val & val >=magx[i+1][j+1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}

	else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	  if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	  if( magx[i-1][j+1]<val & val >=magx[i+1][j-1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	 
      }
    }
  }
  

  //write_f_mat2(smot, "edge_nms");
			
  count = vi.size();
  int *ii, *jj, num, tr,tc;
  ii = new int [count+1];
  jj = new int [count+1];
  
  int R[8]={0,1,0,-1,  1, 1,-1,-1}; // right hand out edge follow
  int C[8]={1,0,-1,0,  1,-1,-1, 1};  // canny edges rarely branch.
  int R2[8]={1,0,-1,0, 1,-1,-1,1},     C2[8]={0,1,0,-1, 1,1,-1,-1};
  //
  for(i=0; i<count; i++){
    
    if( smot[vi[i] ][vj[i] ] > high){ // begin edge follow 8-connected!
      num = 1;
      ii[0] = vi[i]; jj[0] = vj[i];
      smot[ ii[0] ][ jj[0] ] = -4;

      for( int k=0; k<num; k++){  // right hand out edge follow.
	int rx=ii[k], cx=jj[k];
	for( int z=0; z<8; z++){
	  tr = rx+R[z];
	  tc = cx+C[z];
	  if(smot[tr][tc] >0){
	    ii[num] = tr; jj[num] = tc;
	    smot[ tr ][ tc ] = -1;
	    num ++;
	    break;
	  }
	}
      }
      // reverse followed stuff and go off in the other direction.
      Exf_chain eching;
      eching.pts.clear();
      Exf ex;
      for(int q=num-1; q>=0; q--){
	ex.r = ii[q];
	ex.c = jj[q];
	ex.dr = drx[ii[q]][jj[q]];
	ex.dc = dcx[ii[q]][jj[q]];
	eching.pts.push_back(ex);
	//edge[ii[q]][jj[q]] = 1;
      }
      
      //ii[0] = ii[num-1]; jj[0] = jj[num-1];
      num = 1;
      
      for( int k=0; k<num; k++){  // right hand out edge follow, walking backwards.
	int rx=ii[k], cx=jj[k];
	for( int z=0; z<8; z++){
	  tr = rx+R2[z];
	  tc = cx+C2[z];
	  if(smot[tr][tc] >0){
	    ii[num] = tr; jj[num] = tc;
	    smot[ tr ][ tc ] = -2;
	    num ++;
	    break;
	  }
	}
      }
      
      
      if( eching.pts.size() + num > 7){
	for(int q=1; q<num; q++){
	  ex.r = ii[q];
	  ex.c = jj[q];
	  ex.dr = drx[ii[q]][jj[q]];
	  ex.dc = dcx[ii[q]][jj[q]];
	  eching.pts.push_back(ex);
	  //edge[ii[q]][jj[q]] = 1;
	}
	eching.num = chains.size();
	//chains.push_back(eching);
        
	num = eching.pts.size();
	for(int q=0; q<num; q++){
	  tr = eching.pts[q].r;
	  tc = eching.pts[q].c;
	  edge[tr][tc] = 1;
	}

	// turning angle.
	num = eching.pts.size();
	float tangles[num];
	float ang = atan2( eching.pts[0].dr, eching.pts[0].dc), ang2, dang;
	theta = ang;
	i = eching.pts[0].r;
	j = eching.pts[0].c;
	if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	  float y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r;
	  eching.pts[0].mc = (float) eching.pts[0].c + x;
	}
	else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mc = (float) eching.pts[0].c + x;

	}

	else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	 float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  eching.pts[0].mc = (float) eching.pts[0].c;
	}
	else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	  float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mc = (float) eching.pts[0].c + x;
  	}



	for(int q=1; q<num; q++){
	  i = eching.pts[q].r;
	  j = eching.pts[q].c;
	  ang2 = atan2( eching.pts[q].dr, eching.pts[q].dc);
	  dang = ang2-ang;
	
	  if( dang >=  6.283185307)
	    dang = -6.283185307 + dang;
	  if( dang <=  -6.283185307)
	    dang = 6.283185307 + dang;

	  if( dang >= 3.141592654 )
	    dang = -6.283185307 + dang;
	  if( dang <= -3.141592654 )
	    dang = 6.283185307 +dang;

	  tangle[i][j] = dang;
	  tangles[q] = dang;
	  ang = ang2;
	  eching.pts[q].cvtr = dang;

	  theta = ang;
	  if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	    float y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r;
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }
	  else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }
	  
	  else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    eching.pts[q].mc = (float) eching.pts[q].c;
	  }
	  else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }

	}
	if(  (fabs(eching.pts[num-1].r - eching.pts[0].r) <=2 ) && (fabs(eching.pts[num-1].c - eching.pts[0].c) <=2 )){
	  ang = atan2( eching.pts[0].dr, eching.pts[0].dc);
	  dang = ang-ang2;
	
	  if( dang >=  6.283185307)
	    dang = -6.283185307 + dang;
	  if( dang <=  -6.283185307)
	    dang = 6.283185307 + dang;

	  if( dang >= 3.141592654 )
	    dang = -6.283185307 + dang;
	  if( dang <= -3.141592654 )
	    dang = 6.283185307 +dang;

	  tr = eching.pts[0].r;
	  tc = eching.pts[0].c;
	  tangles[0] = dang;
	  tangle[tr][tc] = dang;

	  eching.pts[0].cvtr = dang;

	}
	else{ tangles[0] = 0;}


	// find high curvature curve regions and lines.
	//	HC_pt_line_split( eching, tangles, num, hc_pts, lines); 
	std::vector<Curve> curves;
	Exf_vector_split( eching.pts, lines, hc_pts, curves);

	chains.push_back(eching);

      }
    }
  }

  r = magx[0];
  pend = r + nr*nc;
  for( ; r<pend; r++ ){
    *r = sqrt( *r );
  }
  //write_f_mat2(smot, "edge_fol");
  // write_f_mat2(tangle, "tangle");
   
  // edge.write_rlm_ascii("edge.rlm");
	
  delete [] ii;
  delete [] jj;
  //	delete [] strength_hist;
  
  return(count);
}

/*
  break up chain into high curvature bits and lines


*/

void HC_pt_line_split( Ex_chain &chain, float *tangles, int num, std::vector<Edgel> &hc_pts, std::vector<Line> &lines)
{
  
  int id=0, id2=0, on = 0;
  float sum, n;
  Edgel pt;
  std::vector<PT> cnrs;
  int rs[num], cs[num];

   PT idxx(0, 0 );
   cnrs.push_back(idxx);


  for(int k=0; k<num; k++){
    rs[k]=chain.pts[k].r;
    cs[k]=chain.pts[k].c;
    if( on == 0){
      if( fabs(tangles[k]) > 0.15 ){  // change of state from line to curve
	id = k;
	on = 1;
      }
    }
    else{  
      if( fabs(tangles[k]) < 0.15 ){ // change of state  from corner to line
	on = 0;
	n = (float) k-id;
	float mr=0, mc=0, dr=0, dc=0;
	sum = 0;
	for( int z=id; z<k; z++){    // corner location.
	  sum += fabs(tangles[z]);
	  mr += (float) chain.pts[z].r;
	  mc += (float) chain.pts[z].c;

	  dr += chain.pts[z].dr;
	  dc += chain.pts[z].dc;
	}
	if( sum > 0.3 ){
	  id2 = k;                   // end of corner pixel patch.
	  mr /= n;	mc /= n;
	  dr /= n;	dc /= n;
	  pt.mr = mr;	pt.mc = mc;
	  pt.tr = dr;	pt.tc = dc;
	  pt.cvtr = sum;
	  hc_pts.push_back(pt);

	  PT idx(id,id2);
	  cnrs.push_back(idx);
	}
      }
    }
          
  }
  if( on > 0 ){
    int k=num-1;
    n = (float) k-id;
    float mr=0, mc=0, dr=0, dc=0;
    sum = 0;
    for( int z=id; z<k; z++){    // corner location.
      sum += fabs(tangles[z]);
      mr += (float) chain.pts[z].r;
      mc += (float) chain.pts[z].c;
      
      dr += chain.pts[z].dr;
      dc += chain.pts[z].dc;
    }
    if( sum > 0.3 ){
      id2 = k;                   // end of corner pixel patch.
      mr /= n;	mc /= n;
      dr /= n;	dc /= n;
      pt.mr = mr;	pt.mc = mc;
      pt.tr = dr;	pt.tc = dc;
      pt.cvtr = sum;
      hc_pts.push_back(pt);
      
      PT idx(id,id2);
      cnrs.push_back(idx);
    }
    else{
      PT idx(num-1, 0 );
      cnrs.push_back(idx);
    }
  
  }
  else{
    PT idx(num-1, 0 );
    cnrs.push_back(idx);
  }
  

  int *p1, *p2;
  int  numl = lines.size();

  int numc = cnrs.size();
  if( numc > 0 ){
    numc = cnrs.size();
    for( int q=1; q<numc; q++){
      int id1 = cnrs[q-1].r, id2 = cnrs[q].r;
      p1 = &(rs[id1]);
      p2 = &(cs[id1]);
      int numx = 1+id2-id1;
      if( numx > 3 )
	splitx2( p1, p2, numx, id1, lines);
    }

    int numl2 = lines.size();
    int nums = numl2-numl;
    for(;numl<numl2; numl++){
      lines[numl].id = chain.num;
    }

  }
  
  return;
}

/*
   routine to go from corners to features with scale and orientaiton.

   uses total turning angle in a pseudo octagon as a means of setting scale.


 */

void corner_scale( Matrix<float> &tangle, std::vector<PT> &corns,  std::vector<Corner> &feats,
		   Matrix<float> &dr, Matrix<float> &dc, Matrix<float> &mag)
{
  int nr=tangle.rows(),nc =tangle.cols(), r,c, r0,c0;
  feats.clear();
  
  Matrix<int> zone(nr,nc);
  Matrix<int> tbins(nr,nc);
  // Matrix<float> thetas(nr,nc);

  
  int *pint;
  float *ptr, *pend, *pdr, *pdc;
  ptr = dr[0];
  pdc = dc[0];
  pint = tbins[0];
  pend = ptr + nr*nc;
  //pdr = thetas[0];
  for( ; ptr<pend; ptr++, pdc++, pint++){//, pdr++){
    float theta = atan2( *ptr, *pdc)*18/3.141592654 +17;
    //*pdr = theta;
    *pint = (int) round(theta);
    if( *pint < 0 )
      *pint = 35;

  }

  //  write_mat(thetas, "thetas.mut");

  int nump = corns.size(), odd = 1, k, num, num2;
  int rr[3*(nr+nc)], cc[3*(nr+nc)], rr2[3*(nr+nc)], cc2[3*(nr+nc)];
  int R[8] = {1,0,-1,0, 1,1,-1,-1}, C[8] = {0,1,0,-1, 1,-1,-1,1};
  
  for(int p=0; p<nump; p++){
    int val = p+1;
    float tang = 0;
    rr[0] = corns[p].r;
    cc[0] = corns[p].c;
    Corner crn(corns[p].r,corns[p].c);

    zone[rr[0]][cc[0]] = val;
    num = 1;
    odd = 1;
    for(k=0; k<100; k++){
      if( odd > 0 ){
	num2 = 0;
	for( int q=0; q<num; q++){
	  r0 = rr[q];
	  c0 = cc[q];
	  for( int g=0; g<4; g++){
	    r = r0 + R[g];
	    c = c0 + C[g];
	    if( r <nr && r >=0 && c < nc && c >= 0 ){
	      if( zone[r][c] != val){
		tang += fabs( tangle[r][c]);
		zone[r][c] = val;
		rr2[num2] = r;
		cc2[num2] = c;
		crn.ohist[ tbins[r][c] ]+= mag[r][c];
		num2 ++;
	      }
	    }
	  }
	}
      }
      else{
	num = 0;
	for( int q=0; q<num2; q++){
	  r0 = rr2[q];
	  c0 = cc2[q];
	  for( int g=0; g<8; g++){
	    r = r0 + R[g];
	    c = c0 + C[g];
	    if( r <nr && r >=0 && c < nc && c >= 0 ){
	      if( zone[r][c] != val){
		tang += fabs( tangle[r][c]);
		zone[r][c] = val;
		rr[num] = r;
		cc[num] = c;
		crn.ohist[ tbins[r][c] ]+= mag[r][c];
		num ++;
	      }
	    }
	  }
	}
      }
      //fprintf(stderr,"%d  %d num %d   num2 %d,  tang %.3f \n", p, k, num, num2, tang); 

      if( tang > 12.5664 ){
	break;
      }
      odd *= -1;
    }
    crn.v = (float) k;
    std::vector<float> thetas;
    int numt = find_thetas( crn.ohist, 36, thetas);
    // if( numt > 1){
    //fprintf(stderr," thetas %d, id %d\n", numt, p );
    //}

    for( int th=0; th<numt; th++){
      crn.theta = thetas[th];
      feats.push_back( crn );
    }
    corns[p].v = (float) k;
  }

  //  write_mat(zone, "zone.mat");

  return;
}

/* find the dominant direction in the nhbd of the corner feature.
   

 */
int find_thetas(float *h, int num, std::vector<float> &thetas)
{
  float pix = 3.141592654/18.0;
  float thx[36] = { -3.141592654, -17*pix, -16*pix, -15*pix, -14*pix, -13*pix, -12*pix, -10*pix, -9*pix, 
		    -8*pix, -7*pix, -6*pix, -5*pix, -4*pix, -3*pix, -2*pix, -pix, 0.0, 
		    pix, 2*pix, 3*pix, 4*pix, 5*pix, 6*pix, 7*pix, 8*pix, 9*pix, 
		    10*pix, 11*pix, 12*pix, 13*pix, 14*pix, 15*pix, 16*pix, 17*pix};

  float sum = 0.0, mx = 0.0;
  int id=-1;
  PT hx;
  std::vector<PT> hxs;
  for (int i=0; i<num; i++){
    sum += h[i];
    hx.v = h[i];
    hx.r = i;
    hxs.push_back( hx );
  }
  thetas.clear();
  sort(hxs.begin(), hxs.end());
  thetas.push_back( thx[hxs[0].r]);
  num = 1;

  for (int i=1; i<5; i++){
    if( hxs[i].v > sum/12.0){   // > 3 * the mean.
      thetas.push_back( thx[hxs[i].r]);
      num++;
    }
  }

  return(num);
}


void write_corner_ohist( std::vector<Corner> &m, char * fname )
{
   int nr = m.size(), nc = 4 + 36;
   FILE *fp = fopen(fname, "w");
    
   fprintf(fp,"%d %d\n",nr, nc);

   for (int i=0; i<nr; i++){
     fprintf(fp,"%d %d %.1f %.2f", m[i].r, m[i].c, m[i].v, m[i].theta );
      for (int j=0; j<36; j++){
	  fprintf(fp," %.1f", m[i].ohist[j]);
      }
      fprintf(fp,"\n");
   }

   fclose(fp);
   
   return;
}


void write_corner_ohist_groh( std::vector<Corner> &m, char * fname )
{
   int nr = m.size(), nc = 4 + 36 + 17*16;
   FILE *fp = fopen(fname, "w");
    
   fprintf(fp,"%d %d\n",nr, nc);

   for (int i=0; i<nr; i++){
     fprintf(fp,"%d %d %.1f %.2f", m[i].r, m[i].c, m[i].v, m[i].theta );
      for (int j=0; j<36; j++){
	  fprintf(fp," %.1f", m[i].ohist[j]);
      }
      fprintf(fp,"   ");
      for (int j=0; j<17*16; j++){
	  fprintf(fp," %.1f", m[i].groh[j]);
      }

      fprintf(fp,"\n");
   }

   fclose(fp);
   
   return;
}


/* 
   find native area type features about the corner locations.

   done by resampling the original image then taking derivatives rahter than 
   pyramidal rescaling

   this may or may not be the right thing to do.
 */


void area_feature( Matrix<float> &im, std::vector<Corner> &feats)
{
  int nump = feats.size(), nr=im.rows(), nc=im.cols();

  Matrix<float> pat(39,39);
  
  // make the rescaleable sampling grid, do this as a bilinear thing?
  int m[39][39] = {
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,10,10,10,10,10,10,10,10,10,10,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,9,9,10,10,10,10,10,10,10,10,10,10,10,11,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,9,9,9,9,9,9,9,9,2,2,2,2,2,2,2,2,2,11,11,11,11,11,11,11,11,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,9,9,9,9,9,9,9,1,1,1,2,2,2,2,2,2,2,3,3,3,11,11,11,11,11,11,11,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,9,9,9,9,9,9,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,11,11,11,11,11,11,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,9,9,9,9,9,9,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,11,11,11,11,11,11,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,9,9,9,9,9,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,11,11,11,11,11,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,9,9,9,1,1,1,1,1,1,0,0,0,0,0,0,0,3,3,3,3,3,3,11,11,11,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,1,1,1,1,1,0,0,0,0,0,0,0,0,0,3,3,3,3,3,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,1,1,0,0,0,0,0,0,0,0,0,0,0,3,3,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,8,8,7,7,0,0,0,0,0,0,0,0,0,0,0,5,5,4,4,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,16,16,8,7,7,7,7,7,0,0,0,0,0,0,0,0,0,5,5,5,5,5,4,12,12,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,16,16,16,15,15,15,7,7,7,7,7,7,0,0,0,0,0,0,0,5,5,5,5,5,5,13,13,13,12,12,12,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,15,15,15,15,15,7,7,7,7,7,7,7,6,6,6,6,6,5,5,5,5,5,5,5,13,13,13,13,13,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,15,15,15,15,15,15,7,7,7,7,7,7,6,6,6,6,6,5,5,5,5,5,5,13,13,13,13,13,13,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,15,15,15,15,15,15,7,7,7,7,6,6,6,6,6,6,6,5,5,5,5,13,13,13,13,13,13,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,15,15,15,15,15,15,15,7,7,7,6,6,6,6,6,6,6,5,5,5,13,13,13,13,13,13,13,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,15,15,15,15,15,15,15,15,6,6,6,6,6,6,6,6,6,13,13,13,13,13,13,13,13,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,13,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,13,13,13,13,13,13,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,15,15,15,15,14,14,14,14,14,14,14,14,14,14,14,13,13,13,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,15,15,14,14,14,14,14,14,14,14,14,14,14,13,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,14,14,14,14,14,14,14,14,14,14,14,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};
  
  float er,ec, pr,pc, r,c, r0, c0, rp, rm, cp,cm;

  for(int p=0; p<nump; p++){
    pat.clear();
    Corner &corn = feats[p];
    float scal = corn.v/ 15.0;
    r0 = corn.r;
    c0 = corn.c;

    // theta represents the orientation of the feature patch.
    ec = cos(corn.theta) * scal;  
    er = sin(corn.theta) * scal;   // should be the normal to the curve....
    pr = -ec; pc = er;

    for(int R=-16; R<17; R++){
      for(int C=-16; C<17; C++){
	r = r0 + R*er + C*pr;
	c = c0 + R*ec + C*pc;     // local  bilinear sample.
	rp = ceilf(r);	rm = floorf(r);
	cp = ceilf(c);	cm = floorf(c);

	if( rp > 0 && cp > 0 && rp<nr && cp < nc){
	  if(( rp == r ) && (cp != c )){
	    pat[R+19][C+19] = fabs(c-cm)*im[ (int) round(r)][(int) cp]  + fabs(c-cp)*im[ (int) round(r)][(int) cm] ;
	  }
	  if(( rp != r ) && (cp == c )){
	    pat[R+19][C+19] = fabs(r-rm)*im[(int) rp][(int) round(c)]  + fabs(r-rp)*im[(int) rm][(int) round(c)] ;
	  }
	  if(( rp != r ) && (cp != c )){
	    pat[R+19][C+19] = fabs(r-rm)*fabs(c-cm)*im[ (int) rp][(int) cp] + fabs(r-rp)*fabs(c-cm)*im[ (int) rm][(int) cp] + 
	      fabs(r-rm)*fabs(c-cp)*im[ (int) rp][(int) cm] + fabs(r-rp)*fabs(c-cp)*im[ (int) rm][(int) cm];
	  }
	  if(( rp == r ) && (cp == c )){
	    pat[R+19][C+19] = im[ (int) round(r) ][(int) round(c)];
	  }
	}
      }
    }
    //  now process and sample the patch.
    /*  if( p < 20 ){
	char nam[128];
	sprintf(nam,"pat%d.mut",p);
	write_mat( pat, nam );
	} */

    for(int R=4; R<35; R++){
      for(int C=4; C<35; C++){
	if( m[R][C] >= 0 ){
	  float dr = pat[R+1][C] - pat[R-1][C];
	  float dc = pat[R][C+1] - pat[R][C-1];
	  float mag = sqrt(dr*dr + dc*dc);
	  float thet = atan2(dr,dc);
	  int tbin = (int) round(thet *8/3.141592654 + 7);
	  if( tbin < 0 ) tbin = 15;
	  corn.groh[tbin + m[R][C]*16] += mag;
	}
      }
    } 
  }

  return;
}

/* split chian into line segments

 */

void splitx2( int *i, int *j, int num, int offset, std::vector<Line> &lines)
{
  int x=0,z=0, g = 0;
  int si=0;
  int sjj=0;
  int idx=0;
  float max=0.0, dist= 0.0;
  Line line;
  int *splif;

  splif = new int [ num+10 ];
  memset((void *) splif, 0, (num + 10) * sizeof(int));
  splif[num-1] = 1;

  // find most distant point from start
  si = i[0]; 
  sjj = j[0];
  idx = num-1;
  for( int k=1; k<num; k++){
    dist = (i[k]-si)*(i[k]-si) + (j[k]-sjj)*(j[k]-sjj);

    if( dist > max ){
      max = dist;
      idx = k;
    }
  }
  splif[num-1] = 1;
  int start=0, stop=0;

  start = 0;
  stop = idx;
  splif[stop] = 1;


curse:;
  float h = peak_dist2(i, j, start, stop, &idx);
  if( h > 9.0 ){
    stop = idx;
    splif[stop] = 1;

    goto curse;
  }
  else {
    if( stop - start >= 3 ){
      line.r1 = i[start];
      line.c1 = j[start];
      line.r2 = i[stop];
      line.c2 = j[stop];
      line.v1 = start+offset;
      line.v2 = stop+offset;
      
      /*
      line.s = (float) i[(int) ((start+stop)/2 )];
      line.w = (float) j[(int) ((start+stop)/2 )];
      
      //     line.pc =   line.r2 - line.r1; // perpendiculat to line.
      //    line.pr = - line.c2 + line.c1;
      line.pc =   line.c2 - line.c1; // parallel to line.
      line.pr =   line.r2 + line.r1;
      float mag = sqrt( line.pc*line.pc + line.pr*line.pr);
      line.length = mag;
      if( mag > 0.0 ){
         line.pr /= mag;
         line.pc /= mag;
	 } */
      lines.push_back( line );
    }
    start = stop;
    // look along splif for next 1
    for( int k=start+1; k<num; k++){
       if( splif[k] > 0 ){
          stop = k;
          goto curse;
       }
    }
  }

  delete [] splif;
  return;

}


float peak_dist2(int *i, int *j, int start, int stop, int *idx)
{

  // find the most distant point from the line (si,sj) -> (ei, ej)

  float h,dr,dc, mag;
  float si, sj, tr, tc;
  int id;
  float max = 0.0;

  si = (float) i[start];
  sj = (float) j[start];
  
  dr = i[stop] - si;
  dc = j[stop] - sj;
  mag = sqrt( dr*dr + dc*dc );
  if( mag > 0 ){
    dr/=mag; 
    dc/=mag;
    for( int k=start; k<=stop; k++){
      tr = (i[k]-si); tc = (j[k]-sj);
      mag = tr*tr + tc*tc;
      float tmp = dr*tr + dc*tc;
      h = mag - tmp*tmp;
      if( h> max ){
	max = h;
	id = k;
      }
    }
  }
  *idx = id;
  return( max );

}


float dmax_Exf( std::vector<Exf> &exs, int x1, int x2, int &x3)
{
  float mr = exs[x1].mr, mc = exs[x1].mc;
  float dr = exs[x2].mr - mr,  dc = exs[x2].mc - mc;
  float dmax = -1;
  
  float mag = sqrt(dr*dr + dc*dc);
  if( mag > 0 ){
    dr = dr/mag;
    dc = dc/mag;

    for(int k=x1+1; k<x2; k++){
      float h =  fabs((exs[k].mr - mr)*dc  - (exs[k].mc - mc)*dr);
      if( h > dmax ){
	dmax = h;
	x3 = k;
      }
    }
  }
  else{
    for(int k=x1+1; k<x2; k++){
      float h =  sqrt((exs[k].mr - mr)*(exs[k].mr - mr)  - (exs[k].mc - mc)*(exs[k].mc - mc));
      if( h > dmax ){
	dmax = h;
	x3 = k;
      }
    }
  }

  return(dmax);
}

/* 
 * break up an edgel chain that has orientation information present.
 *
 * note this normalises the [dr,dc] to mag = 1.
 *
 */

void Exf_chain_split( std::vector<Exf> &pts, int start, int stop, std::vector<Line> &lines, int chain_id )
{
   int num = stop-start+1, nums=0;
   float dr,dc,mag, r,c, tr,tc;
   std::vector<int> ids;
   //int *rs = new int [num];
   //int *cs = new int [num];
   float rs[num], cs[num];
   float  drs[num], dcs[num], cvtrs[num];
   
   for(int k=0; k<num; k++){
     int q=start+k;
     rs[k] = pts[q].mr;
     cs[k] = pts[q].mc;
     mag = sqrt( pts[q].dr* pts[q].dr + pts[q].dc* pts[q].dc);
     drs[k] = ((float) pts[q].dr)/mag;
     dcs[k] = ((float) pts[q].dc)/mag;
     //pts[q].dr = drs[k];
     //pts[q].dc = dcs[k];
   }
   
   // compute cvtr first. (num > 10)
   cvtrs[0]=0;
   for(int k=1; k<num; k++){
     cvtrs[k]=1-fabs( drs[k]* drs[k-1] + dcs[k]*dcs[k-1]);
   }
   
   ids.push_back(0);
   for(int k=2; k<num-2; k++){
     if( cvtrs[k] > 0.0341 ){ 
       float f = cvtrs[k];
       if(( f > cvtrs[k-1] ) &&  (f > cvtrs[k-2] ) &&  (f >= cvtrs[k+1] ) &&  (f >= cvtrs[k+2] ) ){  
	 ids.push_back(k);
       }
     }
   }
   ids.push_back(num-1);
   nums = ids.size();
   // split up sub chains into lines!
   
   float *p1, *p2;
   int  numl = lines.size();
   for(int k=0; k<nums-1; k++){
     int id1 = ids[k], id2 = ids[k+1];
     p1 = &(rs[id1]);
     p2 = &(cs[id1]);
     int numx = 1+id2-id1;
     splitx( p1, p2, numx, id1, lines);
     
   }
   int numl2 = lines.size();
   nums = numl2-numl;
   for(;numl<numl2; numl++){
     lines[numl].id = chain_id;
   }
   
   //delete rs;
   //delete cs;
   return;
}







void splitx( float *i, float *j, int num, int offset, std::vector<Line> &lines)
{
  int x,z;
  x = 0;
  z = 0;
  int g = 0;

  int si=0;
  int sjj=0;
  int idx=0;
  float max=0.0, dist= 0.0;
  Line line;
  int *splif;

  splif = new int [ num+10 ];
  memset((void *) splif, 0, (num + 10) * sizeof(int));
  splif[num-1] = 1;

  // find most distant point from start
  si = i[0]; 
  sjj = j[0];
  idx = num-1;
  for( int k=1; k<num; k++){
    dist = (i[k]-si)*(i[k]-si) + (j[k]-sjj)*(j[k]-sjj);

    if( dist > max ){
      max = dist;
      idx = k;
    }
  }
  splif[num-1] = 1;
  int start=0, stop=0;

  start = 0;
  stop = idx;
  splif[stop] = 1;


curse:;
  float h = peak_dist(i, j, start, stop, &idx);
  if( h > 9.0 ){
    stop = idx;
    splif[stop] = 1;

    goto curse;
  }
  else {
    if( stop - start >= 3 ){
      line.r1 = i[start];
      line.c1 = j[start];
      line.r2 = i[stop];
      line.c2 = j[stop];
      line.v1 = start+offset;
      line.v2 = stop+offset;
      

      line.s = (float) i[(int) ((start+stop)/2 )];
      line.w = (float) j[(int) ((start+stop)/2 )];
      
      //      line.pc =   line.r2 - line.r1; // perpendiculat to line.
      //    line.pr = - line.c2 + line.c1;
      line.pc =   line.c2 - line.c1; // parallel to line.
      line.pr =   line.r2 + line.r1;
      float mag = sqrt( line.pc*line.pc + line.pr*line.pr);
      line.length = mag;
      if( mag > 0.0 ){
         line.pr /= mag;
         line.pc /= mag;
      }
      lines.push_back( line );
    }
    start = stop;
    // look along splif for next 1
    for( int k=start+1; k<num; k++){
       if( splif[k] > 0 ){
          stop = k;
          goto curse;
       }
    }
  }

  delete [] splif;
  return;

}


float peak_dist(float *i, float *j, int start, int stop, int *idx)
{

  // find the most distant point from the line (si,sj) -> (ei, ej)

  float h,dr,dc, mag;
  float si, sj, tr, tc;
  int id;
  float max = 0.0;

  si = (float) i[start];
  sj = (float) j[start];
  
  dr = i[stop] - si;
  dc = j[stop] - sj;
  mag = sqrt( dr*dr + dc*dc );
  if( mag > 0 ){
    dr/=mag; 
    dc/=mag;
    for( int k=start; k<=stop; k++){
      tr = (i[k]-si); tc = (j[k]-sj);
      mag = tr*tr + tc*tc;
      float tmp = dr*tr + dc*tc;
      h = mag - tmp*tmp;
      if( h> max ){
	max = h;
	id = k;
      }
    }
  }
  *idx = id;
  return( max );

}





/* 
  not a particualrly vintage way of breaking things into lines.

 */
void split_line_Exf( std::vector<Exf> &exs, int x1, int x2, std::vector<Line> &lines)
{
  int flag = 1, end = x2, x3;
  float thresh;
  if( x2-x1 < 30 ){
    thresh = 1 + 60/(x2-x1);
  }
  else{
    thresh = 3;
  }
  
  while( flag > 0){
    float h = dmax_Exf( exs, x1, x2, x3);
    
    if( h > thresh ){
      x2 = x3;
    }
    else{
      Line lin;
      lin.r1 = exs[x1].mr;
      lin.r2 = exs[x2].mr;
      lin.c1 = exs[x1].mc;
      lin.c2 = exs[x2].mc;
      lines.push_back(lin);
      if( x2 == end ){
	flag = 0;
      }
      else{
	x1 = x2;
	x2 = end;
      }
    }
  }
  return;
}




void Exf_vector_split( std::vector<Exf> &exs, std::vector<Line> &lines, std::vector<Exf> &pts, std::vector<Curve> &curves){

  int nump = exs.size(), id=0, id1;
  float mag;

  for(int p=1; p<nump-1; p++){
    if( (fabs(exs[p].cvtr) >  fabs(exs[p-1].cvtr)) &&  (fabs(exs[p].cvtr) >=  fabs(exs[p+1].cvtr)) ){
      if( fabs(exs[p].cvtr) +  fabs(exs[p-1].cvtr) + fabs(exs[p+1].cvtr) > 0.5236 ){
	id1 = p-1;
	if( id1-id >= 3){
	  split_line_Exf( exs, id, id1, lines);
	}
	id = p+1;

	float xr=0.0,xc=0.0;
	for(int k=id1; k<=id; k++){
	  xr += exs[k].mr;
	  xc += exs[k].mc;
	}
	Exf hcp;
	hcp.mr = xr/3;
	hcp.mc = xc/3;
	hcp.r = id;
	hcp.c = id1;
	pts.push_back( hcp);
      }
    }
  }
 
  if( nump-id > 4 ){
     split_line_Exf( exs, id, nump-1, lines);
  }


  return;
}



/* start from grey image to keep things simple
 * one scale edge detection works well enough.
 *
 *
 * returns chains and edge map and turning angle map.
 * 
 * edges are given to sub pixel accuracy (quadratic fit) to reduce alyassing woes.
 *
 * image is assumed to be lightly smoothen before start.

 */

int canny_grey(Matrix<float> &im, float low, float high, BWImage &edge,
	       Matrix<float> &drx, Matrix<float> &dcx, std::vector<Exf> &edgels )
{
  int i,j, *ptr, *ptc;
  float mag, *pend, *r, *ptf, dr,dc;
  edgels.clear();
		
  int nr = im.rows(); int br = nr-1;
  int nc = im.cols(); int bc = nc-1;
		
  Matrix<float> magx(nr,nc);
  /* Matrix<float> drx(nr,nc);
     Matrix<float> dcx(nr,nc); */
  Matrix<float> smot(nr,nc);
  

  // compute derivatives
  for(i=1; i<br; i++){
    r = im[i];
    pend = r + bc;
    r += 1;
    
    for( j=1; r<pend; r++, j++ ){
      dcx[i][j] = (*(r+1 )) - (*(r-1 ));
      drx[i][j] = (*(r+nc)) - (*(r-nc));
      magx[i][j] = dcx[i][j]*dcx[i][j] + drx[i][j]*drx[i][j];
    }
  }
  //write_f_mat2(magx, "edge_mag");
  //write_f_mat2(drx, "drx");
  //write_f_mat2(dcx, "dcx");
		
  edge.clear();
  smot.clear();
	
  // non-maximum supression
  int count;
  float val, theta;
  //  float p1 = -7.0*3.141592654/8.0, p2=-5.0*3.141592654/8.0, p3=-3*3.141592654/8.0, p4=-3.141592654/8.0, p5=3.141592654/8.0, p6=3.0*3.141592654/8.0, p7=5.0*3.141592654/8.0, p8=7.0*3.141592654/8.0;
  float p1 = -6.8*3.141592654/8.0, p2=-5.2*3.141592654/8.0, p3=-2.8*3.141592654/8.0, p4=-1.2*3.141592654/8.0, 
        p5 =  1.2*3.141592654/8.0, p6=2.8*3.141592654/8.0,  p7= 5.2*3.141592654/8.0, p8= 6.8*3.141592654/8.0;
  float nms_strength=0, vr, vc, x, r2 = 1.4142;
  vector<int> vi, vj;
	
  count = 0;
  for(i=1; i<br; i++){
    ptf = magx[i];
    pend  = ptf + bc;
    ptf  += 1;
    
    for(int j=1 ; ptf<pend; ptf++, j++){
      if( *ptf > low ){
	vr = (float) drx[i][j];
	vc = (float) dcx[i][j];
	val = *ptf;
	
	theta = atan2(vr,vc);
	if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	  if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	  if( magx[i-1][j-1]<val & val >=magx[i+1][j+1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}

	else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	  if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	  if( magx[i-1][j+1]<val & val >=magx[i+1][j-1] ){
	    smot[i][j] = val;
	    vi.push_back(i);
	    vj.push_back(j);
	  }
	}
	 
      }
    }
  }
  

  //write_f_mat2(smot, "edge_nms");
			
  count = vi.size();
  int *ii, *jj, num, tr,tc, cnt =0;
  ii = new int [count+1];
  jj = new int [count+1];
  

  //  fprintf(stderr,"of edge things %d\n", count);

  int R[8]={0,1,0,-1,  1, 1,-1,-1}; // right hand out edge follow
  int C[8]={1,0,-1,0,  1,-1,-1, 1};  // canny edges rarely branch.
  int R2[8]={1,0,-1,0, 1,-1,-1,1},     C2[8]={0,1,0,-1, 1,1,-1,-1};
  //
  for(i=0; i<count; i++){
    
    if( smot[vi[i] ][vj[i] ] > high){ // begin edge follow 8-connected!
      num = 1;
      ii[0] = vi[i]; jj[0] = vj[i];
      smot[ ii[0] ][ jj[0] ] = -4;

      for( int k=0; k<num; k++){  // right hand out edge follow.
	int rx=ii[k], cx=jj[k];
	for( int z=0; z<8; z++){
	  tr = rx+R[z];
	  tc = cx+C[z];
	  if(smot[tr][tc] >0){
	    ii[num] = tr; jj[num] = tc;
	    smot[ tr ][ tc ] = -1;
	    num ++;
	    break;
	  }
	}
      }
      // reverse followed stuff and go off in the other direction.
      Exf_chain eching;
      eching.pts.clear();
      Exf ex;
      for(int q=num-1; q>=0; q--){
	ex.r = ii[q];
	ex.c = jj[q];
	ex.dr = drx[ii[q]][jj[q]];
	ex.dc = dcx[ii[q]][jj[q]];
	eching.pts.push_back(ex);
	//edge[ii[q]][jj[q]] = 1;
      }
      
      //ii[0] = ii[num-1]; jj[0] = jj[num-1];
      num = 1;
      
      for( int k=0; k<num; k++){  // right hand out edge follow, walking backwards.
	int rx=ii[k], cx=jj[k];
	for( int z=0; z<8; z++){
	  tr = rx+R2[z];
	  tc = cx+C2[z];
	  if(smot[tr][tc] >0){
	    ii[num] = tr; jj[num] = tc;
	    smot[ tr ][ tc ] = -2;
	    num ++;
	    break;
	  }
	}
      }
      
      
      if( eching.pts.size() + num > 7){
	cnt ++;                              // this is a chain index id.
	for(int q=1; q<num; q++){
	  ex.r = ii[q];
	  ex.c = jj[q];
	  ex.dr = drx[ii[q]][jj[q]];
	  ex.dc = dcx[ii[q]][jj[q]];
	  eching.pts.push_back(ex);
	  //edge[ii[q]][jj[q]] = 1;
	}
	eching.num = cnt;
	//chains.push_back(eching);
        
	num = eching.pts.size();
	for(int q=0; q<num; q++){
	  tr = eching.pts[q].r;
	  tc = eching.pts[q].c;
	  edge[tr][tc] = 1;
	}

	// turning angle.
	num = eching.pts.size();
	float tangles[num];
	float ang = atan2( eching.pts[0].dr, eching.pts[0].dc), ang2, dang;
	theta = ang;
	i = eching.pts[0].r;
	j = eching.pts[0].c;
	if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	  float y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r;
	  eching.pts[0].mc = (float) eching.pts[0].c + x;
	}
	else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mc = (float) eching.pts[0].c + x;

	}

	else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	 float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  eching.pts[0].mc = (float) eching.pts[0].c;
	}
	else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	  float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mr = (float) eching.pts[0].r + x;
	  y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	  x = (y-y0)/((y+y0)*2+0.0001);
	  eching.pts[0].mc = (float) eching.pts[0].c + x;
  	}


	for(int q=1; q<num; q++){
	  i = eching.pts[q].r;
	  j = eching.pts[q].c;
	  theta = atan2( eching.pts[q].dr, eching.pts[q].dc);
	 
	  if( (theta <= p1 || theta >= p8) || ( theta >= p4 && theta <= p5)){
	    float y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r;
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }
	  else if( (theta >= p1 && theta <= p2) || ( theta >= p5 && theta <= p6)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }
	  
	  else if( (theta >= p2 && theta <= p3) || ( theta >= p6 && theta <= p7)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    eching.pts[q].mc = (float) eching.pts[q].c;
	  }
	  else if( (theta >= p3 && theta <= p4) || ( theta >= p7 && theta <= p8)){
	    float y = sqrt(magx[i+1][j]), y0 = sqrt(magx[i-1][j]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mr = (float) eching.pts[q].r + x;
	    y = sqrt(magx[i][j+1]), y0 = sqrt(magx[i][j-1]); 
	    x = (y-y0)/((y+y0)*2+0.0001);
	    eching.pts[q].mc = (float) eching.pts[q].c + x;
	  }

	}
	/* make 3 poxel long edgels.
	   as features for the texture model.
	*/

	Exf edg;
	int q, cnt2 = 0;
	edg.r = cnt;
	for( q=2; q<num; q+=3){
	  edg.mr = (eching.pts[q-2].mr +eching.pts[q-1].mr +eching.pts[q].mr)/3.0;
	  edg.mc = (eching.pts[q-2].mc +eching.pts[q-1].mc +eching.pts[q].mc)/3.0;

	  edg.dr = (eching.pts[q-2].dr +eching.pts[q-1].dr +eching.pts[q].dr)/3.0;
	  edg.dc = (eching.pts[q-2].dc +eching.pts[q-1].dc +eching.pts[q].dc)/3.0;

	  edg.c = cnt2;
	  cnt2 ++;
	  
	  edgels.push_back(edg);
	}
	
	//else{ tangles[0] = 0;}
	

	// find high curvature curve regions and lines.
	//	HC_pt_line_split( eching, tangles, num, hc_pts, lines); 
	//	std::vector<Curve> curves;
	//Exf_vector_split( eching.pts, lines, hc_pts, curves);
	
	//chains.push_back(eching);
	
      }
    }
  }
  
  /*  r = magx[0];
      pend = r + nr*nc;
      for( ; r<pend; r++ ){
      *r = sqrt( *r );
      } */
  //write_f_mat2(smot, "edge_fol");
  // write_f_mat2(tangle, "tangle");
   
  // edge.write_rlm_ascii("edge.rlm");
	
  delete [] ii;
  delete [] jj;
  //	delete [] strength_hist;
  
  return(count);
}

/* build two hop feature pairwise relations for histogramming

   note affects(normalises) fx.gradient info.

 */


void texture_feats( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  std::vector<Adj> adjs;
  Adj adj;
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }


  int numf = fx.size();
  float thetas[numf], pie = 18.0/3.141592654;

  for( int k=0; k<numf; k++){
    adj.k = k;
    adjs.push_back(adj);
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;
  }

  int numl = links.size();
  for( int i=0; i<numl; i++){
    Link &x = links[i];
    if( x.a >=0 && x.b >= 0 ){
      adjs[x.a].ids.insert(x.b);
      adjs[x.b].ids.insert(x.a);
    }

    /* if( x.a >=0 && x.c >= 0 ){
      adjs[x.a].ids.insert(x.c);
      adjs[x.c].ids.insert(x.a);
    }

    if( x.c >=0 && x.b >= 0 ){
      adjs[x.b].ids.insert(x.c);
      adjs[x.c].ids.insert(x.b);
      } */
  }
  fprintf(stderr,"adjacencies created \n");

  std::vector<Adj> adjs2(adjs);
  
  for( int i=0; i<numf; i++){
    Adj &aj = adjs[i];
    // fprintf(stderr,"%d %d ",i, aj.ids.size());

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label;
      // add neighbours of neighbours of aj  as neighbours of adjs2[i]
      adjs2[i].ids.insert( adjs[b].ids.begin(), adjs[b].ids.end() ); 
      ++label;

      // fprintf(stderr,"  %d   ", b );
    }

    //fprintf(stderr,"%d %d\n", aj.ids.size(), adjs2[i].ids.size());

  }
  fprintf(stderr,"adjacencies up hopped \n");

  //for( int i=0; i<numf; i++){
  //  fprintf(stderr,"%d  %d  %d\n",i, adjs[i].ids.size(),  adjs2[i].ids.size());
  //}


  // now stuff things into the histograms.
  float dr, dc, dth, mr,mc, br,bc;
  int R,C, bin, tbin;

  for( int i=0; i<numf; i++){
    Adj &aj = adjs2[i];
    Exf &g = fx[i];
    
    mr = g.mr;
    mc = g.mc;
    float numa = (float) aj.ids.size();    
    float numax = 1.0/numa;
    br = (3.0*mr)/nr;
    bc = (3.0*mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    float theta = thetas[i];
    tbin = (int) round( theta) + 18;
    if( tbin >=36 ) tbin = 0;

    hists[bin].ohist[tbin] ++;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label, dbin;
      dr = fx[b].mr - mr;
      dc = fx[b].mc - mc;
      
      float dx = dr*g.dr + dc*g.dc;
      if( dx < -10 )
	dbin =0;
      else 
	if( dx < -5)
	  dbin = 1;
	else
	  if( dx < -1.5)
	    dbin = 2;
	  else
	    if( dx <= 1.5)
	      dbin = 3;
	    else
	      if( dx <= 5)
		dbin = 4;
	      else
	  	if( dx <= 10)
		  dbin = 5;
		else
		  dbin = 6;
		  
      dth = theta - thetas[b];
      if( dth >= 36 ) dth -= 36;
      if( dth <= -36 ) dth += 36;
      
      if( dth > 18 ) dth -= 36;
      if( dth < -18 ) dth += 36;

      int dtbin =  (int) round( dth) + 18;
      if( tbin >=36 ) tbin = 0;

      hists[bin].woh[dbin*36+dtbin] += numax;
      label++;

    }

  }


  for( int h=0; h<9; h++){
    for(int q=0; q<36; q++){
      hists[h].r += hists[h].ohist[q];
    }
  }


  return;
}




/* straight pairwise geometric histogram

   max interaction length 30 pixels.

   3 histograms
   orientation, cross distance vs angle
   parallel distance vs included angle as pairs.

 */


void texture_pgh( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  std::vector<Adj> adjs;
  Adj adj;
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }


  int numf = fx.size();
  float thetas[numf], pie = 18.0/3.141592654, mug;

  for( int k=0; k<numf; k++){
    adj.k = k;
    adjs.push_back(adj);
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;
  }

  int numl = links.size();
  for( int i=0; i<numl; i++){
    Link &x = links[i];
    if( x.a >=0 && x.b >= 0 ){
      mug = (fx[x.a].mr - fx[x.b].mr)* (fx[x.a].mr - fx[x.b].mr) + (fx[x.a].mc - fx[x.b].mc)* (fx[x.a].mc - fx[x.b].mc);
      if( mug < 900){
	adjs[x.a].ids.insert(x.b);
	adjs[x.b].ids.insert(x.a);
      }
    }

  }
  // fprintf(stderr,"adjacencies created \n");

  int taken[numf], lx[numf];
  memset( (void *) taken, 0, numf*sizeof( int ) );

  std::vector<Adj> adjs2(adjs);  
  for( int i=0; i<numf; i++){
    Adj &aj = adjs[i];
    // fprintf(stderr,"%d %d ",i, aj.ids.size());
    int nx=0;
    taken[i] = i;

    float mr = fx[i].mr, mc = fx[i].mc;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label;
      lx[nx] = b;
      taken[b] = i;
      nx++;
      ++label;

      // fprintf(stderr,"  %d   ", b );
    }
    int j=0;
    while( j<nx){
      Adj &aj = adjs[lx[j]];
      // fprintf(stderr,"%d %d ",i, aj.ids.size());
      
      std::set<int> :: const_iterator label=aj.ids.begin();
      while(label != aj.ids.end() ){
	int b = *label;
	if( taken[b] != i ){
	  taken[b] = i;
	  mug = (mr - fx[b].mr)* (mr - fx[b].mr) + (mc - fx[b].mc)* (mc - fx[b].mc);
	  if( mug < 900 ){
	    adjs2[i].ids.insert(b);
	    lx[nx] = b;
	    nx++;
	  }
	}
	++label;
      }
      j++;
    }


    //fprintf(stderr,"%d %d\n", aj.ids.size(), adjs2[i].ids.size());

  }
  //fprintf(stderr,"adjacencies up hopped \n");
  

  //for( int i=0; i<numf; i++){
  //  fprintf(stderr,"%d  %d  %d\n",i, adjs[i].ids.size(),  adjs2[i].ids.size());
  //}


  // now stuff things into the pairwise geometric histograms.
  float dr, dc, dth, mr,mc, br,bc;
  int R,C, bin, tbin;

  for( int i=0; i<numf; i++){
    Adj &aj = adjs2[i];
    Exf &g = fx[i];
    
    mr = g.mr;
    mc = g.mc;
    float numa = (float) aj.ids.size();    
    float numax = 10.0/numa;
    br = (3.0*mr)/nr;
    bc = (3.0*mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    float theta = thetas[i];
    tbin = (int) round( theta) + 18;
    if( tbin >=36 ) tbin = 0;

    hists[bin].ohist[tbin] ++;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label, dbin;
      dr = fx[b].mr - mr;
      dc = fx[b].mc - mc;
      
      float dx = dr*g.dr + dc*g.dc;
      if( dx < -25 )
	dbin =0;
      else 
      if( dx < -20 )
	dbin =1;
      else 
      if( dx < -15 )
	dbin =2;
      else 
      if( dx < -9 )
	dbin =3;
      else 
      if( dx < -5 )
	dbin =4;
      else 
      if( dx < -2)
	dbin =5;
      else 
      if( dx < -1 )
	dbin =6;
      else 
      if( dx <= 1 )
	dbin =7;
      else 
      if( dx <= 2 )
	dbin =8;
      else 
	if( dx <= 5)
	  dbin = 9;
	else
	  if( dx <= 9 )
	    dbin = 10;
	  else
	    if( dx <= 15)
	      dbin = 11;
	    else
	      if( dx <= 20)
		dbin = 12;
	      else
	  	if( dx <= 25)
		  dbin = 13;
		else
		  dbin = 14;
      
      dth = theta - thetas[b];
      if( dth >= 36 ) dth -= 36;
      if( dth <= -36 ) dth += 36;
      
      if( dth > 18 ) dth -= 36;
      if( dth < -18 ) dth += 36;
      
      int dtbin =  (int) round( dth) + 18;
      if( tbin >=36 ) tbin = 0;
      
      hists[bin].woh[dbin*36+dtbin] += numax;
      if( fabs(dx) < 10 ){
	dx = fabs(dc*g.dr - dr*g.dc);
	if( dx < 2 )
	  dbin =0;
	else 
	  if( dx < 5)
	    dbin = 1;
	  else
	    if( dx < 9)
	      dbin = 2;
	    else
	      if( dx < 14)
		dbin = 3;
	      else
	      if( dx < 19)
		dbin = 4;
	      else
	  	if( dx < 24)
		  dbin = 5;
		else
		  dbin = 6;

	float arse = fabs( g.dr * fx[b].dr +  g.dc * fx[b].dc );
	if( arse >= 1.0) arse = 0.999999;
	dth = acos(arse)*(36.0/3.141592636);
	dtbin =  (int) round( dth);
	if( dtbin >= 18 ) dtbin = 17;

	hists[bin].pgh[dbin*18+dtbin] += numax;
      }
      label++;
      
    }

  }


  for( int h=0; h<9; h++){
    for(int q=0; q<36; q++){
      hists[h].v += hists[h].ohist[q];
    }
  }


  return;
}


/* 

// two pass distance transform
// edges >0, distance as -ve integers.
// record the orientation of the nearest edge?


*/
void  odtrans(BWImage &im, BWImage &shad)
{
  int nr = im.rows(), nr2 = nr-2, nc = im.cols(), nc2 = nc-2;
  int v, ms, mx, m=-nr-nc, nc1 = nc+1, nx1 = nc-1;
  
  for (register unsigned int r = 1; r < nr - 1; r++) {
    int *p = im[r] + 1;
    int *ps = shad[r]+1;
    
    for (int *pend = p + nc2; p < pend; p++, ps++)
      if (*p == 0){
         mx = m;
	 ms = -1;

         v = *(p - nc1);
         if (v > 0){ *p= -1; *ps = *(ps-nc1); goto jump1;}
         else if ((v < 0) && (v > mx)) { mx = v;   ms = *(ps-nc1); }
                  
         v = *(p - nc);
         if (v > 0){    *p= -1; *ps = *(ps-nc); goto jump1;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nc);  }
  
         v = *(p - nx1);
         if (v > 0){     *p= -1;  *ps = *(ps-1);  goto jump1;}
         else if ((v < 0) && (v > mx)){    mx = v;  ms = *(ps-1); }

         v = *(p - 1);
         if (v > 0) {   *p= -1; *ps = *(ps-1); goto jump1;}
         v = *(p+1);
         if (v > 0) {   *p= -1; *ps = *(ps+1); goto jump1;}
         v = *(p + nx1);
         if (v > 0){    *p= -1; *ps = *(ps+nx1); goto jump1;}
         v = *(p + nc);
         if (v > 0){    *p= -1; *ps = *(ps+nc); goto jump1;}
         v = *(p + nc1);
         if (v > 0){     *p= -1; *ps = *(ps+nc1); goto jump1;}

         *p = mx-1;
	 *ps = ms;
         jump1:;
      }
  }
  for (register unsigned int j = nc2; j >= 1; --j){
    int *p = im[nr2] + j;
    int *ps = shad[nr2] + j;
    
    for (register int i = nr2; i >= 1; --i, p -= nc, ps -= nc)
      if (*p <= 0){
         mx = m;
	 ms = -1;
         
         v = *(p + nc1);
         if (v > 0){     *p= -1; *ps = *(ps+nc1); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps+nc1); }
         
         v = *(p + nc);
         if (v > 0){    *p= -1; *ps = *(ps+nc); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps+nc); }
                           
         v = *(p + nx1);
         if (v > 0){    *p= -1; *ps = *(ps+nx1); goto jump2;}
         else if ((v < 0) && (v > mx))  {  mx = v; ms = *(ps+nx1); }
         
         v = *(p - 1);
         if (v > 0) {   *p= -1; *ps = *(ps-1); goto jump2;}
         else if (( v < 0) && (v > mx)) {   mx = v; ms = *(ps-1); }
         
         v = *(p - nc);
         if (v > 0){    *p= -1; *ps = *(ps-nc); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nc); }
  
         v = *(p - nx1);
         if (v > 0){     *p= -1; *ps = *(ps-nx1); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nx1); }

         v = *(p+1);
         if (v > 0) {   *p= -1; *ps = *(ps+1); goto jump2;}
         else if ((v < 0) && (v > mx)){    mx = v; ms = *(ps+1); }

         v = *(p - nc1);
         if (v > 0){ *p= -1; *ps = *(ps-nc1); goto jump2;}
         else if ((v < 0) && (v > mx)){  mx = v; ms = *(ps-nc1); }

  
         *p = mx-1;
	 *ps = ms;
         jump2:;
      }
  }
  
  return;
}


/* 
  Make an orientation edge map grom the edge map.

  Here we will use 11 oreintations and 4 edge labels.

*/


void vorx(BWImage &edge, Matrix<float> &drx, Matrix<float> &dcx, BWImage &orient )
{
  int nr = edge.rows(); int br = nr-1;
  int nc = edge.cols(); int bc = nc-1;

  int *p=edge[0], *ps = orient[0];
  float *pr = drx[0], *pc = dcx[0];

  int *pend  = p+nr*nc;
     
  for( ; p<pend; p++, ps++, pr++, pc++){
    if( *p > 0 ){
      float theta = atan2(*pr, *pc);
      if( theta < 0 ) theta += 3.141592636;
      theta *= 12.0/3.141592636;
      int tbin = (int) (round(theta)); 
      if( tbin > 11 ) tbin = 0;
      *ps = tbin;
    }
  }

  // shaddow boundary orientation set.

  ps = orient[0]; 
  pend = ps + nc;
  p = orient[br];
  for(; ps<pend; p++, ps++){
    *ps = 12;
    *p = 13;
  }

  ps = orient[0]; pend = ps + nc*nr;
 
  for(; ps<pend; ps++){
    *ps = 14;
    ps+=bc;
    *ps  = 15;
  }

  // edge boundary instantiation
  ps = edge[0]; 
  pend = ps + nc;
  p = edge[br];
  for(; ps<pend; p++, ps++){
    *ps = 1;
    *p = 1;
  }

  ps = edge[0]; pend = ps + nc*nr;
 
  for(; ps<pend; ps++){
    *ps = 1;
    ps+=bc;
    *ps  = 1;
  }

  return;
}

/* 

// 6 pass distance transform
// edges >0, distance as -ve integers.
// record the orientation of the nearest edge?


*/
void  odtrans6(BWImage &im, BWImage &shad)
{
  int nr = im.rows(), nr2 = nr-2, nc = im.cols(), nc2 = nc-2;
  int v, ms, mx, m=-nr-nc, nc1 = nc+1, nx1 = nc-1;
  
  // 2 horizontal pass to start.

  for ( int r = 1; r < nr - 1; r++) {
    int *p = im[r] + 1;
    int *ps = shad[r]+1;
    int *px = im[r];
    
    for (int *pend = p + nc2; p < pend; p++, ps++){
      if (*p == 0){
	v = *(p - 1);
	if (v > 0){ *p= -1; *ps = *(ps-1); }
	else if (v < 0){ *p = *(p-1)-1; *ps = *(ps-1);}
      }
    }
    // reverse pass
    p--; ps--;
    for (; p > px; p--, ps--){
      if (*p < 0){  // look right.
	if( *(p+1) > 0 ){
	  *p = -1;
	  *ps = *(ps+1);
	}
	else{
	  if (*p < *(p+1) ) { *p = *(p+1)-1; *ps = *(ps+1);}
	}	
      }
    }
    
  }
  
  //shad.write_rlm_ascii("orient1.rlm");
  //im.write_rlm_ascii("dtrans1.rlm");
  
  // 2 vertical passes.
  for (int c=0; c<nc; c++){
    int *p = &(im[1][c]);
    int *ps = &(shad[1][c]);
    for(int r2=1; r2<=nr2; r2++, p+=nc, ps+=nc){
      if( *p < 0 ){  // look up
	if( *(p-nc) > 0){
	  *p = -1; *ps = *(ps-nc);
	}
	else if( *p < *(p-nc) ){
	  *p = *(p-nc) -1; *ps = *(ps-nc);
	}
	
      }
    }

    p = &(im[nr2][c]);
    ps = &(shad[nr2][c]);
    
    for(int r2=nr2; r2>0; r2--, p-=nc, ps-=nc){
      if( *p < 0 ){  // look down
	if( *(p+nc) > 0){
	  *p = -1; *ps = *(ps+nc);
	}
	else if( *p < *(p+nc) ){
	  *p = *(p+nc) -1; *ps = *(ps+nc);
	}
	
      }
    }


  }

  //shad.write_rlm_ascii("orient2.rlm");
  //im.write_rlm_ascii("dtrans2.rlm");


  return;


  for (int r = 1; r < nr - 1; r++) {
    int *p = im[r] + 1;
    int *ps = shad[r]+1;
    
    for (int *pend = p + nc2; p < pend; p++, ps++)
      if (*p == 0){
         mx = m;
	 ms = -1;

         v = *(p - nc1);
         if (v > 0){ *p= -1; *ps = *(ps-nc1); goto jump1;}
         else if ((v < 0) && (v > mx)) { mx = v;   ms = *(ps-nc1); }
                  
         v = *(p - nc);
         if (v > 0){    *p= -1; *ps = *(ps-nc); goto jump1;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nc);  }
  
         v = *(p - nx1);
         if (v > 0){     *p= -1;  *ps = *(ps-1);  goto jump1;}
         else if ((v < 0) && (v > mx)){    mx = v;  ms = *(ps-1); }

         v = *(p - 1);
         if (v > 0) {   *p= -1; *ps = *(ps-1); goto jump1;}
         v = *(p+1);
         if (v > 0) {   *p= -1; *ps = *(ps+1); goto jump1;}
         v = *(p + nx1);
         if (v > 0){    *p= -1; *ps = *(ps+nx1); goto jump1;}
         v = *(p + nc);
         if (v > 0){    *p= -1; *ps = *(ps+nc); goto jump1;}
         v = *(p + nc1);
         if (v > 0){     *p= -1; *ps = *(ps+nc1); goto jump1;}

         *p = mx-1;
	 *ps = ms;
         jump1:;
      }
  }
  for (register unsigned int j = nc2; j >= 1; --j){
    int *p = im[nr2] + j;
    int *ps = shad[nr2] + j;
    
    for (register int i = nr2; i >= 1; --i, p -= nc, ps -= nc)
      if (*p <= 0){
         mx = m;
	 ms = -1;
         
         v = *(p + nc1);
         if (v > 0){     *p= -1; *ps = *(ps+nc1); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps+nc1); }
         
         v = *(p + nc);
         if (v > 0){    *p= -1; *ps = *(ps+nc); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps+nc); }
                           
         v = *(p + nx1);
         if (v > 0){    *p= -1; *ps = *(ps+nx1); goto jump2;}
         else if ((v < 0) && (v > mx))  {  mx = v; ms = *(ps+nx1); }
         
         v = *(p - 1);
         if (v > 0) {   *p= -1; *ps = *(ps-1); goto jump2;}
         else if (( v < 0) && (v > mx)) {   mx = v; ms = *(ps-1); }
         
         v = *(p - nc);
         if (v > 0){    *p= -1; *ps = *(ps-nc); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nc); }
  
         v = *(p - nx1);
         if (v > 0){     *p= -1; *ps = *(ps-nx1); goto jump2;}
         else if ((v < 0) && (v > mx)) {   mx = v; ms = *(ps-nx1); }

         v = *(p+1);
         if (v > 0) {   *p= -1; *ps = *(ps+1); goto jump2;}
         else if ((v < 0) && (v > mx)){    mx = v; ms = *(ps+1); }

         v = *(p - nc1);
         if (v > 0){ *p= -1; *ps = *(ps-nc1); goto jump2;}
         else if ((v < 0) && (v > mx)){  mx = v; ms = *(ps-nc1); }

  
         *p = mx-1;
	 *ps = ms;
         jump2:;
      }
  }
  
  return;
}






/* New area and orientaiton feature.

   maybe it work, maybe it not work.
 */


void fabulist_pgh( BWImage &dtran, BWImage &orient,  std::vector<Xhist>  &hists)
{
  int nr = dtran.rows(), nc = dtran.cols();
  int br1 = nr/3, br2 = (nr*2)/3, bc1 = nr/3, bc2 = (nr*2)/3;
  
  int h1[9] = {0,0,0,br1,br1,br1,br2,br2,br2};
  int h2[9] = {br1,br1,br1,br2,br2,br2, nr,nr,nr};
  
  int h3[9] = {0,bc1,bc2,0,bc1,bc2,0,bc1,bc2};
  int h4[9] = {bc1,bc2,nc,bc1,bc2,nc,bc1,bc2,nc};
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }

  for( int h=0; h<9; h++){

    Xhist &x = hists[h];
    int a=h3[h],b = h4[h];
    for(int r=h1[h]; r<h2[h]; r++){
      for(int c=a; c<b; c++){
	int tbin = orient[r][c];
	int dx = dtran[r][c], dbin;
	if( dx < -100 )
	  dbin =11;
	else 
	  if( dx < -70 )
	    dbin =10;
	  else 
	    if( dx < -50 )
	      dbin =9;
	    else 
	      if( dx < -35 )
		dbin =8;
	      else 
		if( dx < -25 )
		  dbin =7;
		else 
		  if( dx <= -16)
		    dbin =6;
		  else 
		    if( dx <= -10 )
		      dbin =5;
		    else 
		      if( dx <= -6 )
			dbin =4;
		      else 
			if( dx <= -4 )
			  dbin =3;
			else 
			  if( dx <= -2)
			    dbin = 2;
			  else
			    if( dx < 0 )
			      dbin = 1;
			    else
			      dbin = 0;
	x.shad[dbin*16+tbin] ++;
      }
    }
    
  }


  return;
}







/* straight pairwise geometric histogram

   max interaction length 30 pixels.

 */


void texture_pgh2( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  std::vector<Adj> adjs;
  Adj adj;
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }


  int numf = fx.size();
  float thetas[numf], pie = 18.0/3.141592654, mug;

  for( int k=0; k<numf; k++){
    adj.k = k;
    adjs.push_back(adj);
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;
  }

  int numl = links.size();
  for( int i=0; i<numl; i++){
    Link &x = links[i];
    if( x.a >=0 && x.b >= 0 ){
      mug = (fx[x.a].mr - fx[x.b].mr)* (fx[x.a].mr - fx[x.b].mr) + (fx[x.a].mc - fx[x.b].mc)* (fx[x.a].mc - fx[x.b].mc);
      if( mug < 900){
	adjs[x.a].ids.insert(x.b);
	adjs[x.b].ids.insert(x.a);
      }
    }

  }
  // fprintf(stderr,"adjacencies created \n");

  int taken[numf], lx[numf];
  memset( (void *) taken, 0, numf*sizeof( int ) );

  std::vector<Adj> adjs2(adjs);  
  for( int i=0; i<numf; i++){
    Adj &aj = adjs[i];
    // fprintf(stderr,"%d %d ",i, aj.ids.size());
    int nx=0;
    taken[i] = i;

    float mr = fx[i].mr, mc = fx[i].mc;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label;
      lx[nx] = b;
      taken[b] = i;
      nx++;
      ++label;

      // fprintf(stderr,"  %d   ", b );
    }
    int j=0;
    while( j<nx){
      Adj &aj = adjs[lx[j]];
      // fprintf(stderr,"%d %d ",i, aj.ids.size());
      
      std::set<int> :: const_iterator label=aj.ids.begin();
      while(label != aj.ids.end() ){
	int b = *label;
	if( taken[b] != i ){
	  taken[b] = i;
	  mug = (mr - fx[b].mr)* (mr - fx[b].mr) + (mc - fx[b].mc)* (mc - fx[b].mc);
	  if( mug < 900 ){
	    adjs2[i].ids.insert(b);
	    lx[nx] = b;
	    nx++;
	  }
	}
	++label;
      }
      j++;
    }


    //fprintf(stderr,"%d %d\n", aj.ids.size(), adjs2[i].ids.size());

  }
  //fprintf(stderr,"adjacencies up hopped \n");
  

  //for( int i=0; i<numf; i++){
  //  fprintf(stderr,"%d  %d  %d\n",i, adjs[i].ids.size(),  adjs2[i].ids.size());
  //}


  // now stuff things into the pairwise geometric histograms.
  float dr, dc, dth, mr,mc, br,bc;
  int R,C, bin, tbin;

  for( int i=0; i<numf; i++){
    Adj &aj = adjs2[i];
    Exf &g = fx[i];
    
    mr = g.mr;
    mc = g.mc;
    float numa = (float) aj.ids.size();    
    float numax = 10.0/numa;
    br = (3.0*mr)/nr;
    bc = (3.0*mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    float theta = thetas[i];
    tbin = (int) round( theta) + 18;
    if( tbin >=36 ) tbin = 0;

    hists[bin].ohist[tbin] ++;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label, dbin;
      dr = fx[b].mr - mr;
      dc = fx[b].mc - mc;
      
      float dx = dr*g.dr + dc*g.dc;   // across the edgel
      float dp = -dr*g.dc + dc*g.dr;  // parallel to Mr edgel.

      float arse = fabs( g.dr * fx[b].dr +  g.dc * fx[b].dc );
      if( arse >= 1.0) arse = 0.999999;
      dth = acos(arse)*(36.0/3.141592636);
      int dtbin =  (int) round( dth);
      if( dtbin >= 18 ) dtbin = 17;
      
      if( fabs(dp) <= 10 ){

	  if( dx < -25 )
	    dbin =0;
	  else 
	    if( dx < -20 )
	      dbin =1;
	    else 
	      if( dx < -15 )
		dbin =2;
	      else 
		if( dx < -9 )
		  dbin =3;
		else 
		  if( dx < -5 )
		    dbin =4;
		  else 
		    if( dx < -2)
		      dbin =5;
		    else 
		      if( dx < -1 )
			dbin =6;
		      else 
			if( dx <= 1 )
			  dbin =7;
			else 
			  if( dx <= 2 )
			    dbin =8;
			  else 
			    if( dx <= 5)
			      dbin = 9;
			    else
			      if( dx <= 9 )
				dbin = 10;
			      else
				if( dx <= 15)
				  dbin = 11;
				else
				  if( dx <= 20)
				    dbin = 12;
				  else
				    if( dx <= 25)
				      dbin = 13;
				    else
				      dbin = 14;
            
	  hists[bin].woh[dbin*18+dtbin] += numax;
	}
	if( fabs(dx) <= 10 ){
	  float dpx = fabs(dp);
	  if( dpx < 2 )
	    dbin =0;
	  else 
	    if( dpx < 5)
	      dbin = 1;
	    else
	      if( dpx < 9)
		dbin = 2;
	      else
		if( dpx < 14)
		  dbin = 3;
		else
		  if( dpx < 19)
		    dbin = 4;
		  else
		    if( dpx < 24)
		      dbin = 5;
		    else
		      dbin = 6;
	  hists[bin].pgh[dbin*18+dtbin] += numax;
	}
	
	// remaining 12 outer quadrants bins to be filled with stuff.
	if( (dx > 10) && (dx <= 20 ) && (dp >10) && (dp <= 20)){
	   hists[bin].pgh2[dtbin] += numax;
	}
	else if( (dx < -10) && (dx >= -20 ) && (dp >10) && (dp <= 20)){
	   hists[bin].pgh2[18+dtbin] += numax;
	}
	else if( (dx < -10) && (dx >= -20 ) && (dp < -10) && (dp >= -20)){
	   hists[bin].pgh2[36+dtbin] += numax;
	}
	else if( (dx > 10) && (dx <= 20 ) && (dp < -10) && (dp >= -20)){
	   hists[bin].pgh2[54+dtbin] += numax;
	}

	else if( (dx > 20) && (dp > 10) ){
	   hists[bin].pgh2[72+dtbin] += numax;
	}
	else if( (dx < -20) && (dp > 10) ){
	   hists[bin].pgh2[90+dtbin] += numax;
	}
	else if( (dx < -20) && (dp < -10) ){
	   hists[bin].pgh2[108+dtbin] += numax;
	}
	else if( (dx > 20) && (dp < -10) ){
	   hists[bin].pgh2[126+dtbin] += numax;
	}

	else if( (dx > 10) && (dx <= 20 ) && (dp > 20)){
	   hists[bin].pgh2[144+dtbin] += numax;
	}
	else if( (dx < -10) && (dx >= -20 ) && (dp > 20) ){
	   hists[bin].pgh2[162+dtbin] += numax;
	}
	else if( (dx < -10) && (dx >= -20 ) && (dp < - 20) ){
	   hists[bin].pgh2[180+dtbin] += numax;
	}
	else if( (dx > 10) && (dx <= 20 ) && (dp < - 20)  ){
	   hists[bin].pgh2[198+dtbin] += numax;
	}


	label++;
      
    }

  }

    /*
  for( int h=0; h<9; h++){
    for(int q=0; q<36; q++){
      hists[h].v += hists[h].ohist[q];
    }
    } */


  return;
}


















/* straight pairwise geometric histogram

   max interaction length 30 pixels.

 */


void texture_pgh3( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  std::vector<Adj> adjs;
  Adj adj;
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }


  int numf = fx.size();
  float thetas[numf], pie = 36.0/3.141592654, pie2 = 36.0/3.141592654;

  for( int k=0; k<numf; k++){
    adj.k = k;
    adjs.push_back(adj);
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    if( thetas[k] > pie2){
      thetas[k] -= pie;
    }
    if( thetas[k] < -pie2){
      thetas[k] += pie;
    }


    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;
  }

  int numl = links.size();
  for( int i=0; i<numl; i++){
    Link &x = links[i];
    if( x.a >=0 && x.b >= 0 ){
      float mug = (fx[x.a].mr - fx[x.b].mr)* (fx[x.a].mr - fx[x.b].mr) + (fx[x.a].mc - fx[x.b].mc)* (fx[x.a].mc - fx[x.b].mc);
      if( mug < 900){
	adjs[x.a].ids.insert(x.b);
	adjs[x.b].ids.insert(x.a);
      }
    }

  }
  // fprintf(stderr,"adjacencies created \n");

  int taken[numf], lx[numf];
  memset( (void *) taken, 0, numf*sizeof( int ) );

  std::vector<Adj> adjs2(adjs);  
  for( int i=0; i<numf; i++){
    Adj &aj = adjs[i];
    // fprintf(stderr,"%d %d ",i, aj.ids.size());
    int nx=0;
    taken[i] = i;

    float mr = fx[i].mr, mc = fx[i].mc;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label;
      lx[nx] = b;
      taken[b] = i;
      nx++;
      ++label;

      // fprintf(stderr,"  %d   ", b );
    }
    int j=0;
    while( j<nx){
      Adj &aj = adjs[lx[j]];
      // fprintf(stderr,"%d %d ",i, aj.ids.size());
      
      std::set<int> :: const_iterator label=aj.ids.begin();
      while(label != aj.ids.end() ){
	int b = *label;
	if( taken[b] != i ){
	  taken[b] = i;
	  float mug = (mr - fx[b].mr)* (mr - fx[b].mr) + (mc - fx[b].mc)* (mc - fx[b].mc);
	  if( mug < 900 ){
	    adjs2[i].ids.insert(b);
	    lx[nx] = b;
	    nx++;
	  }
	}
	++label;
      }
      j++;
    }


    //fprintf(stderr,"%d %d\n", aj.ids.size(), adjs2[i].ids.size());

  }
  //fprintf(stderr,"adjacencies up hopped \n");
  

  //for( int i=0; i<numf; i++){
  //  fprintf(stderr,"%d  %d  %d\n",i, adjs[i].ids.size(),  adjs2[i].ids.size());
  //}


  // now stuff things into the pairwise geometric histograms.
  float dr, dc, dth, mr,mc, br,bc;
  int R,C, bin, tbin;

  for( int i=0; i<numf; i++){
    Adj &aj = adjs2[i];
    Exf &g = fx[i];
    
    mr = g.mr;
    mc = g.mc;
    float numa = (float) aj.ids.size();    
    float numax = 10.0/numa;
    br = (3.0*mr)/nr;
    bc = (3.0*mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    float theta = thetas[i];
    tbin = (int) round( theta) + 18;
    if( tbin >=36 ) tbin = 0;

    //    hists[bin].ohist[tbin] ++;

    std::set<int> :: const_iterator label=aj.ids.begin();
    while(label != aj.ids.end() ){
      int b = *label, dbin;
      dr = fx[b].mr - mr;
      dc = fx[b].mc - mc;
      
      float dx = dr*g.dr + dc*g.dc;   // across the edgel
      float dp = -dr*g.dc + dc*g.dr;  // parallel to Mr edgel.

      float arse = thetas[b] - theta;
      if( arse > pie2){
	arse -= pie;
      }
      if( arse < -pie2){
	arse += pie;
      }
      int dtbin = (int) round( arse + 18 );

      if( dtbin >= 36 ) dtbin = 0;
      if( dtbin < 0 ) dtbin = 36;
      
      if( fabs(dx) <= 10 ){

	if( dp < -24 )
	  dbin =0;
	else 
	  if( dp < -19 )
	    dbin =1;
	  else 
	    if( dp < -14 )
	      dbin =2;
	    else 
	      if( dp < -9 )
		dbin =3;
	      else 
		if( dp < -5 )
		  dbin =4;
		else 
		  if( dp < -2)
		    dbin =5;
		  else 
		    if( dp <= 2 )
		      dbin =6;
		    else 
		      if( dp <= 5)
			dbin = 7;
		      else
			if( dp <= 9 )
			  dbin = 8;
			else
			  if( dp <= 14)
			    dbin = 9;
			  else
			    if( dp <= 19)
			      dbin = 10;
			    else
			      if( dp <= 24)
				dbin = 11;
			      else
				dbin = 12;
	
	  hists[bin].woh[dbin*36+dtbin] += numax;
      }
      if( fabs(dx) > 10 && fabs(dx) <= 20 ){  // outer flank. (partial 180 deg ambiguity).
	
	if( dp < -19 )
	  dbin =0;
	else 
	  if( dp < -14 )
	    dbin =1;
	  else 
	    if( dp < -9 )
	      dbin =2;
	    else 
	      if( dp < -5 )
		dbin =3;
	      else 
		if( dp < -2)
		  dbin =4;
		else 
		  if( dp <= 2 )
		    dbin =5;
		  else 
		    if( dp <= 5)
		      dbin = 6;
		    else
		      if( dp <= 9 )
			  dbin = 7;
		      else
			if( dp <= 14)
			  dbin = 8;
			else
			  if( dp <= 19)
			    dbin = 9;
			  else
			    dbin = 10;

	dtbin = (int) round( arse/2.0 + 9.0 );
	if( dtbin >= 18 ) dtbin = 0;
	if( dtbin < 0 ) dtbin = 17;
	if( dx > 0 )
	  hists[bin].pgh[dbin*18+dtbin] += numax;
	else
	  hists[bin].pgh2[dbin*18+dtbin] += numax;
      }
	

      label++;
      
    }

  }

    /*
  for( int h=0; h<9; h++){
    for(int q=0; q<36; q++){
      hists[h].v += hists[h].ohist[q];
    }
    } */


  return;
}




void Exf_lines( std::vector<Exf> &edg, std::vector<Line> &lines )
{
  int num = edg.size();
  int stop=1, start = 0, val;
  int k=0;

  val = edg[1].r;    // chain id number.
  while( stop<num){
    if( edg[stop-1].r == edg[stop].r ){
      stop ++;
    }
    else{
      if( stop - start >= 5 ){
	Exf_chain_split( edg, start, stop-1, lines, edg[start].r );
      }
      start = stop;
      stop++;
    }
  }
  if( stop - start >= 5 ){
    Exf_chain_split( edg, start, stop-1, lines, edg[start].r );
  }

  return;
}

/* 
   make a line length vs orientation hist set

*/ 


void lines_pgh( std::vector<Xhist> &hists, int nr, int nc, std::vector<Line> &lines )
{
  int num = lines.size();
  float arse = sqrt(nr*nr + nc*nc);  // length bins.
  float lnd = log(arse)/12;

  float nr3 = 3.0/((float) nr);
  float nc3 = 3.0/((float) nc);


  // 10 degree bins.
  float pie = 18.0/3.141592654;


  for( int i=0; i< num; i++){
    Line &lx = lines[i];
    float dr=lx.r2-lx.r1, dc=lx.c2-lx.c1;

    float theta = atan2( dr, dc )*pie; // 180 degree ambiguity
    if( theta < 0){
      theta += 18;
    }
    int dtbin = (int) round( theta);
    
    if( dtbin >= 18 ) dtbin = 0;
    float mag = sqrt(dr*dr + dc*dc);
    int xbin = (int) ( log(mag)/lnd)-3;   // bins linked to image size please.
    if( xbin < 0 ) xbin = 0;
    if( xbin >= 8) xbin = 8;

    
    dr /= 32; dc /= 32;
    float mg = mag/arse;       // size independent length.
    for(int k=0; k<=32; k++){
      float r = (lx.r1+dr*k)*nr3;
      float c = (lx.c1+dc*k)*nc3;
      int C = (int) c;
      int R = (int) r;
      hists[ R*3 + C].pgh2[xbin*18 + dtbin] += mg;  // linear function of line length...
    }

  }
  return;
} 




void texture_ohist( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
  }


  int numf = fx.size();
  float pie = 18.0/3.141592654;


  // now stuff things into the pairwise geometric histograms.
  float dr, dc, dth, mr,mc, br,bc;
  int R,C, bin, tbin;

  for( int i=0; i<numf; i++){
    br = (3.0*fx[i].mr)/nr;
    bc = (3.0*fx[i].mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    float theta =  atan2( fx[i].dr, fx[i].dc)*pie;
    tbin = (int) round( theta) + 18;
    if( tbin >=36 ) tbin = 0;

    hists[bin].ohist[tbin] ++;
  }
  return;
}

/* 
   revised texture model giving cvtr in line and medial axis width off line

*/


void texture_cvtrw( std::vector<Exf> &fx, std::vector<Link> &links, int nr, int nc, std::vector<Xhist>  &hists)
{
  std::vector<Adj> adjs;
  Adj adj;
  /* hists.clear();
  Xhist xh;
  for( int h=0; h<9; h++){
    hists.push_back(xh);
    }*/


  int numf = fx.size(), tbin, bin, R,C;
  float thetas[numf], pie = 18.0/3.141592654, mug, br,bc,
    pie180 = 180.0/3.141592654;

  for( int k=0; k<numf; k++){
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;

    tbin = (int) round( thetas[k]) + 18;
    if( tbin >=36 ) tbin = 0;

    br = (3.0*fx[k].mr)/nr;
    bc = (3.0*fx[k].mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    hists[bin].ohist[tbin] ++;
  }

  float ln15 = log(1.5);
  int numl = links.size();
  for( int i=0; i<numl; i++){
    Link &x = links[i];
    int a = x.a, b=x.b;
    if( a >=0 && b >= 0 ){
      float dr = fx[b].mr - fx[a].mr, dc = fx[b].mc - fx[a].mc;
      mug = dr*dr + dc*dc;
      float mr = (fx[a].mr + fx[b].mr)/2;
      float mc = (fx[a].mc + fx[b].mc)/2;

      float e1 =  fx[a].dr *fx[b].dr + fx[a].dc *fx[b].dc;
      float e2 =  -fx[a].dc *fx[b].dr + fx[a].dr *fx[b].dc;
      float thet = atan2(e2,e1)*pie180;
      float slum = -dr*fx[a].dc + dc*fx[a].dr;
      float slum2 = dr*fx[a].dr + dc*fx[a].dc;

      br = (3.0*mr)/nr;
      bc = (3.0*mc)/nc;
      R = (int) br;
      C = (int) bc;
      bin = R*3+C;
    
      if( fx[a].r == fx[b].r){
	if( abs( fx[a].c - fx[b].c) < 2 ){
	  // genuine +/- curvature of edge with sign of gradient.
	  float swt = slum*thet;

	  //fprintf(stderr, "%.2f ", thet);
	  if( fabs(thet) < 1.0 ) tbin = 0;          // 1 deg
	  else if( fabs(thet) < 2.0 ) 
	    if( swt > 0 )
	      tbin = 1;
	    else
	      tbin = 2;    
	  else if( fabs(thet) < 3.0 ) 
	    if( swt > 0 )
	      tbin = 3;
	    else
	      tbin = 4;    
	  else if( fabs(thet) < 4.0 ) 
	    if( swt > 0 )
	      tbin = 5;
	    else
	      tbin = 6;    
	  else  if( fabs(thet) < 10.0 ) 
	    if( swt > 0 )
	      tbin = 7;
	    else
	      tbin = 8;    
	  else  if( fabs(thet) < 20.0 ) 
	    if( swt > 0 )
	      tbin = 9;
	    else
	      tbin = 10;    
	  else  if( fabs(thet) < 35.0 ) 
	    if( swt > 0 )
	      tbin = 11;
	    else
	      tbin = 12;
	  else  if( fabs(thet) < 65.0 ) 
	    if( swt > 0 )
	      tbin = 13;
	    else
	      tbin = 14;        
	  else if( fabs(thet) < 70.0 ) 
	    if( swt > 0 )
	      tbin = 15;
	    else
	      tbin = 16;  
	  else  if( fabs(thet) < 80.0 ) 
	    if( swt > 0 )
	      tbin = 17;
	    else
	      tbin = 18;    
	  //if( fabs(thet) < 85.0 ) 
	  else
	    if( swt > 0 )
	      tbin = 19;
	    else
	      tbin = 20;
	  hists[bin].woh[13+tbin] ++;  
	}
	else{
	  int wbin;
	  if( mug < 3.0) wbin = 0;
	  else if ( mug < 4.5) wbin = 1;
	  else if ( mug < 6.75) wbin = 2;
	  else if ( mug < 10.125) wbin = 3;
	  else if ( mug < 15.2) wbin = 4;
	  else if ( mug < 22.8) wbin = 5;
	  else if ( mug < 34.2) wbin = 6;
	  else if ( mug < 76.9) wbin = 7;
	  else if ( mug < 115.3) wbin = 8;
	  else if ( mug < 173.0) wbin = 9;
	  else if ( mug < 259.5) wbin = 10;
	  else if ( mug < 389.2) wbin = 11;
	  else wbin = 12;
	  
	  hists[bin].woh[wbin] ++;  // the 13 bins of Aramoth.
	}
      }
      else{
	int wbin;
	if( mug < 3.0) wbin = 0;
	else if ( mug < 4.5) wbin = 1;
	else if ( mug < 6.75) wbin = 2;
	else if ( mug < 10.125) wbin = 3;
	else if ( mug < 15.2) wbin = 4;
	else if ( mug < 22.8) wbin = 5;
	else if ( mug < 34.2) wbin = 6;
	else if ( mug < 76.9) wbin = 7;
	else if ( mug < 115.3) wbin = 8;
	else if ( mug < 173.0) wbin = 9;
	else if ( mug < 259.5) wbin = 10;
	else if ( mug < 389.2) wbin = 11;
	else wbin = 12;
	
	hists[bin].woh[wbin] ++;  // the 13 bins of Aramoth.
      }
    }

  }

  return;
}

/* 
   fit ellipses to the edgels underlying lines.
   
   look for the longest chain of lines segments that can be fitted by an ellispe?
   
   wish to know the total turning angle of the ellipse.
   
*/


void Exf_line_ellipse_fit(std::vector<Exf> &pts, std::vector<Line> &lines, std::vector<Ellipse> &ellipses, int chain_id)
{
  Ellipse ellap;
   
   int numl = lines.size();
   
   int nump = pts.size();
   float rs[nump], cs[nump], drs[nump], dcs[nump];
   float tr[10], tc[10], mr[10], mc[10], e1[5], dr[10], dc[10], mags[5], mdrc;
   int ex;
   
   if( numl >= 4){
      dr[0] = 0; dr[4] = 0; dc[0] = 0; dc[4] = 0;
      for( int k=0; k<nump; k++){
         rs[k] = (float)pts[k].r;
         cs[k] = (float)pts[k].c;
         drs[k] = (float)pts[k].dr;
         dcs[k] = (float)pts[k].dc;
         mdrc = sqrt(drs[k]*drs[k] + dcs[k]*dcs[k] );
         drs[k] /= mdrc;
         dcs[k] /= mdrc;
      }

      for( int x=0; x<=numl-4; x++){
         int id1=lines[x].v1, id2=lines[x+2].v2;
         
         // check line segments for turning angle direction
         for(int e=0; e<4; e++){
            ex = x+e;
            tr[e] = (lines[ex].r2-lines[ex].r1);
            tc[e] = (lines[ex].c2-lines[ex].c1);
            mr[e] = lines[ex].r1;
            mc[e] = lines[ex].c1;
            mags[e] = sqrt( tr[e]*tr[e] + tc[e]*tc[e] );
         }
         if( (mags[0]+mags[1] + mags[2] + mags[3]) > 20){
            
            mr[4] = lines[x+3].r2;         mc[4] = lines[x+3].c2;
            float mug;
            dr[1] = tr[0]/mags[0] + tr[1]/mags[1];
            dc[1] = tc[0]/mags[0] + tc[1]/mags[1];
            mug = sqrt(dr[1]*dr[1] + dc[1]*dc[1]);
            dr[1] /= mug; dc[1] /= mug;
            
            dr[2] = tr[2]/mags[2] + tr[1]/mags[1];
            dc[2] = tc[2]/mags[2] + tc[1]/mags[1];
            mug = sqrt(dr[2]*dr[2] + dc[2]*dc[2]);
            dr[2] /= mug; dc[2] /= mug;
            
            dr[3] = tr[2]/mags[2] + tr[3]/mags[3];
            dc[3] = tc[2]/mags[2] + tc[3]/mags[3];
            mug = sqrt(dr[3]*dr[3] + dc[3]*dc[3]);
            dr[3] /= mug; dc[3] /= mug;
            
            
            
            for(int e=0; e<3; e++){
               e1[e] = (tr[e]*tc[e+1] - tc[e]*tr[e+1])/(mags[e]*mags[e+1]);
            }
            if(( e1[0] > 0 && e1[1] > 0 && e1[2] > 0 ) ||
            (e1[0] < 0 && e1[1] < 0 && e1[2] < 0) ){
               float tx = acos(fabs(e1[0])) + acos(fabs(e1[1])) + acos(fabs(e1[2]));
               if( tx > 1.50 ){
                 
                 float resid;// = ellipse_fit_pseudo9(mr, mc, dr, dc, 5,&(rs[id1]), &(cs[id1]), &(drs[id1]), &(dcs[id1]), id2-id1, ellap);
                  if( ellap.res < 0.2 && ellap.res1 < 6 && resid < 0.3){
                     //std::cerr << "el_fit_err " << ellap.res << " " << ellap.res1 << " " << resid <<  endl;
                     
                     // test for ellipseness
                     if( (4*ellap.A*ellap.C - ellap.B*ellap.B) > 0 ){
                        ellap.flag = chain_id;
                        ellap.idx1 = id1;
                        ellap.idx2 = id2-id1;
                        ellipses.push_back(ellap);
                     }
		     /*                     else{// hyperbola.
                        curv.k = chain_id;
                        curv.id1 = x;
                        curv.dt = acos(e1[0]) + acos(e1[1]) +acos(e1[2]);
                        curv.length = 1+id2-id1;
                        curves.push_back( curv);
			}*/
                  }
               }
            }
         }
      }
   }
	return;
}


/* simple orientaiton 0-18 (with reflection) bins
   and width from edge to edge
   edge widths will have to be normailsed by area....

*/

void texture_OW( std::vector<Exf> &fx, int nr, int nc, std::vector<Xhist>  &hists, BWImage &edge)
{
  std::vector<Adj> adjs;
  Adj adj;

  int numf = fx.size(), tbin, bin, R,C;
  int d, wbin,   wbr1 = nr/10, wbr2 = nr/4, wbr3 = (2*nr)/3;
  int   wbc1 = nc/10, wbc2 = nc/4, wbc3 = (2*nc)/3;
  float thetas[numf], pie = 18.0/3.141592654, mug, br,bc,
    pie180 = 180.0/3.141592654;

  for( int k=0; k<numf; k++){
    thetas[k] = atan2( fx[k].dr, fx[k].dc)*pie;
    float mag = sqrt( fx[k].dr*fx[k].dr + fx[k].dc*fx[k].dc);
    fx[k].dr/= mag;
    fx[k].dc/= mag;

    tbin = (int) round( thetas[k]) + 18;
    if( tbin >=36 ) tbin = 0;

    br = (3.0*fx[k].mr)/nr;
    bc = (3.0*fx[k].mc)/nc;
    R = (int) br;
    C = (int) bc;
    bin = R*3+C;
    
    hists[bin].ohist[tbin] ++;
  }
  
  return;


  int nr1= nr/3, nr2 = (2*nr)/3;
  int nc1= nc/3, nc2 = (2*nc)/3;
  int mc;  

  // width histograms
  for( int r=0; r<nr1; r++){
    int cx = 0;
    
    for( int c=1; c<nc; c++){
      if( edge[r][c] > 0 ){
	mc = (cx+c)/2;
	d = c-cx;
	cx = c;
	bin = 2;
	if( mc < nc1) bin = 0;
	else if( mc <nc2) bin = 1;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin] += d;

      }
    }
    mc = (cx+nc)/2;
    d = nc-cx;
    bin = 2;
    if( mc <nc2) bin = 1;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin] += d; 
  }


  for( int r=nr1; r<nr2; r++){
    int cx = 0;
    
    for( int c=1; c<nc; c++){
      if( edge[r][c] > 0 ){
	mc = (cx+c)/2;
	d = c-cx;
	cx = c;
	bin = 5;
	if( mc < nc1) bin = 3;
	else if( mc <nc2) bin = 4;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin] += d;

      }
    }
    mc = (cx+nc)/2;
    d = nc-cx;
    bin = 5;
    if( mc <nc2) bin = 4;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin] += d; 
  }
  
  for( int r=nr2; r<nr; r++){
    int cx = 0;
    
    for( int c=1; c<nc; c++){
      if( edge[r][c] > 0 ){
	mc = (cx+c)/2;
	d = c-cx;
	cx = c;
	bin = 8;
	if( mc < nc1) bin = 6;
	else if( mc <nc2) bin = 7;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin] += d;

      }
    }
    mc = (cx+nc)/2;
    d = nc-cx;
    bin = 8;
    if( mc <nc2) bin = 7;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin] += d; 
  }


  // vertial strip length histogram.

  // width histograms

  for( int c=0; c<nc1; c++){
    int rx = 0;
    
    for( int r=0; r<nr; r++){
      if( edge[r][c] > 0 ){
	mc = (rx+r)/2;
	d = r-rx;
	rx = r;
	bin = 6;
	if( mc < nr1) bin = 0;
	else if( mc <nr2) bin = 3;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin+7] += d;

      }
    }
    mc = (rx+nr)/2;
    d = nr-rx;
    bin = 6;
    if( mc <nr2) bin = 3;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin+7] += d; 
  }



  for( int c=nc1; c<nc2; c++){
    int rx = 0;
    
    for( int r=0; r<nr; r++){
      if( edge[r][c] > 0 ){
	mc = (rx+r)/2;
	d = r-rx;
	rx = r;
	bin = 7;
	if( mc < nr1) bin = 1;
	else if( mc <nr2) bin = 4;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin+7] += d;

      }
    }
    mc = (rx+nr)/2;
    d = nr-rx;
    bin = 7;
    if( mc <nr2) bin = 4;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin+7] += d; 
  }


  for( int c=nc2; c<nc; c++){
    int rx = 0;
    
    for( int r=0; r<nr; r++){
      if( edge[r][c] > 0 ){
	mc = (rx+r)/2;
	d = r-rx;
	rx = r;
	bin = 8;
	if( mc < nr1) bin = 2;
	else if( mc <nr2) bin = 5;
 
	wbin = 6;
	if( d <=3 ) wbin = 0;
	else if( d <= 6) wbin = 1;
	else if( d <= 12) wbin = 2;
	else if( d <= wbr1) wbin = 3;
	else if( d <= wbr2) wbin = 4;
	else if( d <= wbr3) wbin = 5;

	hists[bin].woh[wbin+7] += d;

      }
    }
    mc = (rx+nr)/2;
    d = nr-rx;
    bin = 8;
    if( mc <nr2) bin = 5;
    wbin = 6;
    if( d <=3 ) wbin = 0;
    else if( d <= 6) wbin = 1;
    else if( d <= 12) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    hists[bin].woh[wbin+7] += d; 
  }



  return;
}




void line_oh( std::vector<Line> &lx, int nr, int nc, std::vector<Xhist>  &hists)
{

  int numl = lx.size();
  float pie = 3.14159265;
  float d = (float) nr;
  if( nc > nr) d = (float) nc;

  float wbr1 = d/10, wbr2 = d/4, wbr3 = (2*d)/3;

  for( int k=0; k<numl; k++){
    float dr = lx[k].r1 -lx[k].r2;
    float dc = lx[k].c1 -lx[k].c2;
    

    float theta = atan2(dr, dc);
    float d = sqrt( dr*dr + dc*dc);
    
    if( theta < 0 ) theta = pie+theta;
    theta = round( theta*4/pie);
    int tbin = (int) theta;
    if( tbin >=4 ) tbin = 0;


    float br = 1.5*(lx[k].r1+ lx[k].r2)/((float)nr);
    float bc = 1.5*(lx[k].c1+ lx[k].c2)/((float)nc);
    int R = (int) br;
    int C = (int) bc;
    int bin = R*3+C;

    int  wbin = 6;
    if( d <=5 ) wbin = 0;
    else if( d <= 10) wbin = 1;
    else if( d <= 20) wbin = 2;
    else if( d <= wbr1) wbin = 3;
    else if( d <= wbr2) wbin = 4;
    else if( d <= wbr3) wbin = 5;
    
    
    hists[bin].pgh[tbin*7+wbin] +=d;
  }

  return;
}

