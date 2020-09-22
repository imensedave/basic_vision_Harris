#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <set>
#include <string>
#include <cfloat>
#include "Imgx.h"
#include "bmp.h"

#include<stdlib.h>




#define CREATOR ""
#define RGB_COMPONENT_COLOR 255

int readPPMx(char *filename, Matrix<short> &grau)//CxImage &cim)
{
  char buff[16];
  //         CxImage img(1,1);
  FILE *fp;
  int c, rgb_comp_color;
  //open PPM file for reading
  fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    return(-1);
  }
  
  //read image format
  if (!fgets(buff, sizeof(buff), fp)) {
    perror(filename);
    fclose(fp);
    return(-1);
  }

  //check the image format
  if (buff[0] != 'P' || buff[1] != '6') {
    fprintf(stderr, "Invalid image format (must be 'P6')\n");
    fclose(fp);
    return(-1);
  }

  // may be #'d comments:(
  int nr=1,nc=1;


  
  //read image size information
  if (sscanf(buff+2, " %d %d %d", &nc, &nr, &rgb_comp_color) != 3) {
    fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
    fclose(fp);
    return(-1);
  }
  
  //check rgb component depth
  if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
    fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
    fclose(fp);
    return(-1);
  }


  //memory allocation for pixel data
  //img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));
  unsigned char *data3 = new unsigned char [nr*nc*3];

  //read pixel data from file
  int n = fread(data3, sizeof(unsigned char), 3 * nr* nc, fp);
  if (n < nr*nc*3) {
    fprintf(stderr, "Error loading image '%s'\n", filename);
    delete [] data3;
    fclose(fp);
    return(-1);
  }

  grau.reshape(nr,nc);

  short *pg = grau[0];
  short *pend = pg+nr*nc;
  unsigned char *pim = data3;

  //  pim++;
  //pim++;

  // seems to be encoded BRG
  
  for( ; pg<pend; pg++){
    int g = (*pim) *4;
    pim ++;
    g += (*pim) *18;
    pim ++;
    g += (*pim) *11;
    pim ++;

    g = g >>5; 
    
    *pg = (short) g;  

    
    //   *pg = *pim;
    // pim+=3;
    

  }

  fclose(fp);
  delete [] data3;
  return(1);
  
}
/*
void writePPM(const char *filename, PPMImage *img)
{
    FILE *fp;
    //open file for output
    fp = fopen(filename, "wb");
    if (!fp) {
         fprintf(stderr, "Unable to open file '%s'\n", filename);
         return(-1);
    }

    //write the header file
    //image format
    fprintf(fp, "P6\n");

    //comments
    fprintf(fp, "# Created by %s\n",CREATOR);

    //image size
    fprintf(fp, "%d %d\n",img->x,img->y);

    // rgb component depth
    fprintf(fp, "%d\n",RGB_COMPONENT_COLOR);

    // pixel data
    fwrite(img->data, 3 * img->x, img->y, fp);
    fclose(fp);
}

void changeColorPPM(PPMImage *img)
{
    int i;
    if(img){

         for(i=0;i<img->x*img->y;i++){
              img->data[i].red=RGB_COMPONENT_COLOR-img->data[i].red;
              img->data[i].green=RGB_COMPONENT_COLOR-img->data[i].green;
              img->data[i].blue=RGB_COMPONENT_COLOR-img->data[i].blue;
         }
    }
    }*/

/*
int main(){
    PPMImage *image;
    image = readPPM("can_bottom.ppm");
    changeColorPPM(image);
    writePPM("can_bottom2.ppm",image);
    printf("Press any key...");
    getchar();
}
 */

int readppm2(char *fname, Matrix<short> &grau)
{
  char s0[513];
    
  std::string line;
  std::string p6("P6");
  std::string octo("#");

  
  std::ifstream myfile;
  myfile.open(fname);

  
  if (myfile.is_open()){

    int nr = 0, nc = 0, bog=0;
   
    getline (myfile,line);
    
    // check string for the string "points"
    int n = (int) line.find( p6);
    if( n == 0 ){
      if( line.length() <= 512){
	copy( line.begin(), line.end(), s0);
	s0[line.length()] = 0;
	int v = sscanf( s0+2, " %d %d %d", &nc, &nr, &bog);
	if( v == 3 ){
	  cerr << "hazah!" << endl;
	}
      
	if( myfile.good() ){
	  int do_next = 1;
	  while ( do_next ){
	    
	    getline (myfile,line);
	    n = (int) line.find( octo); // ### comment line...
	    if( n==0 ){
	      
	    }
	    else{
	      do_next = 0;
	      if( line.length() <= 512){
		copy( line.begin(), line.end(), s0);
		s0[line.length()] = 0;
		int v = sscanf( s0, "%d %d", &nr, &nc);
		if( v == 2 ){
		  
		  cerr << "hazah!" << endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }


  return(1);
}
