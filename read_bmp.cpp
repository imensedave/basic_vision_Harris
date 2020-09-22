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



unsigned char *LoadBitmapFile(char *filename, BITMAPINFOHEADER2 *bitmapInfoHeader)
{
  FILE *filePtr; //our file pointer
  BITMAPFILEHEADER2 bitmapFileHeader; //our bitmap file header
  unsigned char *bitmapImage;  //store image data
  int imageIdx=0;  //image index counter
  unsigned char tempRGB;  //our swap variable
  
  //open filename in read binary mode
  filePtr = fopen(filename,"rb");
  if (filePtr == NULL)
    return NULL;
  
  BMPX bmpx;
  BMPY bmpy;

  fread(&bmpx, sizeof(BMPX),1,filePtr);
  fread(&bmpy, sizeof(BMPY),1,filePtr);


  //read the bitmap file header
  //  fread(&bitmapFileHeader, sizeof(BITMAPFILEHEADER),1,filePtr);
  
  //verify that this is a bmp file by check bitmap id
  if (bmpx.magic[0] != 'B' && bmpx.magic[1] != 'M')//bitmapFileHeader.bfType !=0x4D42)
    {
      fclose(filePtr);
      return NULL;
    }
  
  //read the bitmap info header
  fread(bitmapInfoHeader, sizeof(BITMAPINFOHEADER2),1,filePtr);
  
  //move file point to the begging of bitmap data
  //  fseek(filePtr, bitmapFileHeader.bfOffBits, SEEK_SET);
  fseek(filePtr, bmpy.bmp_offset, SEEK_SET);
  
  //allocate enough memory for the bitmap image data
  int area = (bitmapInfoHeader->biWidth)*(bitmapInfoHeader->biHeight);
  bitmapImage = (unsigned char*)malloc(area);
  if( area == bitmapInfoHeader->biSizeImage || bitmapInfoHeader->biBitCount == 8 ){

  
    //verify memory allocation
    if (!bitmapImage)
      {
	free(bitmapImage);
	fclose(filePtr);
	return NULL;
      }
  
    //read in the bitmap image data
    fread(bitmapImage,area,1,filePtr);
  
    //make sure bitmap image data was read
    if (bitmapImage == NULL)
      {
	fclose(filePtr);
	return NULL;
      }
    
  }

  else if( bitmapInfoHeader->biBitCount == 24 ){
    unsigned char *tmp = (unsigned char*) malloc( area*3);
    
    //verify memory allocation
    if (!bitmapImage)
      {
	free(bitmapImage);
	fclose(filePtr);
	return NULL;
      }
  
    //read in the bitmap image data
    fread(tmp,area*3,1,filePtr);
  
    //make sure bitmap image data was read
    if (bitmapImage == NULL)
      {
	delete []  tmp;
	fclose(filePtr);
	return NULL;
      }

    unsigned char *p= bitmapImage, *pend = p+area, *pd=tmp;
    for(; p<pend; p++){
      pd++;
      *p=*pd;
      pd++;
      pd++;
    }
    delete []  tmp;
  }
  else{
    fclose(filePtr);
    delete [] bitmapImage;
    cerr << "bitmap format not supported " << endl;
    exit(2);
  }





  /* swap the r and b values to get RGB (bitmap is BGR)
     for (imageIdx = 0,imageIdx < bitmapInfoHeader->biSizeImage;imageIdx+=3)
     {
     tempRGB = bitmapImage[imageIdx];
     bitmapImage[imageIdx] = bitmapImage[imageIdx + 2];
     bitmapImage[imageIdx + 2] = tempRGB;
     }
  */
  //close file and return bitmap iamge data
  fclose(filePtr);
  return bitmapImage;
}

