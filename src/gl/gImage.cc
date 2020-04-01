#include "gImage.hh"
#include <stdio.h>
#include <string.h>
GImage::GImage(unsigned int width, unsigned int height):
width(width),
height(height)
{
  // data.resize(width*height*3);
  data = (unsigned char *)malloc(width*height*3);
  memset(data,0,3*width*height);
}

void GImage::setPixel(unsigned int x, unsigned int y, std::array<unsigned int, 3> color)
{
  if(x>=width)return;
  if(y>=height)return;
  unsigned int offset = 3*(y*width+x);
  data[offset+0]=static_cast<unsigned char>(color[2]);
  data[offset+1]=static_cast<unsigned char>(color[1]);
  data[offset+2]=static_cast<unsigned char>(color[0]);
}

void GImage::setPixel(unsigned int x, unsigned int y, unsigned char* color)
{
  if(x>=width)return;
  if(y>=height)return;
  unsigned int offset = 3*(y*width+x);
  data[offset+0]=static_cast<unsigned char>(color[2]);
  data[offset+1]=static_cast<unsigned char>(color[1]);
  data[offset+2]=static_cast<unsigned char>(color[0]);
}

std::array<unsigned int, 3> GImage::getPixel(unsigned int x, unsigned int y) const
{
  std::array<unsigned int, 3> color{0,0,0};
  if(x>=width)return color;
  if(y>=height)return color;
  unsigned int offset = 3*(y*width+x);
  color[0]=data[offset+0];
  color[1]=data[offset+1];
  color[2]=data[offset+2];
  return color;
}

unsigned int GImage::getWidth() const
{
  return width;
}

unsigned int GImage::getHeight() const
{
  return height;
}

void GImage::save(std::string filename) const
{
  unsigned int size =54+ 3*width*height;
unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
unsigned char bmppad[3] = {0,0,0};

bmpfileheader[ 2] = (unsigned char)(size    );
bmpfileheader[ 3] = (unsigned char)(size>> 8);
bmpfileheader[ 4] = (unsigned char)(size>>16);
bmpfileheader[ 5] = (unsigned char)(size>>24);

bmpinfoheader[ 4] = (unsigned char)(       width    );
bmpinfoheader[ 5] = (unsigned char)(       width>> 8);
bmpinfoheader[ 6] = (unsigned char)(       width>>16);
bmpinfoheader[ 7] = (unsigned char)(       width>>24);
bmpinfoheader[ 8] = (unsigned char)(       height    );
bmpinfoheader[ 9] = (unsigned char)(       height>> 8);
bmpinfoheader[10] = (unsigned char)(       height>>16);
bmpinfoheader[11] = (unsigned char)(       height>>24);
FILE *f;
f = fopen(filename.c_str(),"wb");
  // std::ofstream os(filename);
  // if(!os) return;
fwrite(bmpfileheader,1,14,f);
fwrite(bmpinfoheader,1,40,f);

for(unsigned int i=0; i<height; i++)
{
    fwrite(data+(width*(height-i-1)*3),3,width,f);
    fwrite(bmppad,1,(4-(width*3)%4)%4,f);
}

  // os<<"BM"<<std::endl<<width<<" "<<height<<" 255"<<std::endl;
  // for(unsigned int i=0;i<size;i++)
  // {
  //   os<<static_cast<unsigned int>(data[i]);
  //   if( (i+1)%20) os<<' ';
  //   else os<<std::endl;
  // }
  // os<<std::endl;
  fclose(f);
}

