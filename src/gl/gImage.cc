#include "gImage.hh"

GImage::GImage(unsigned int width, unsigned int height):
width(width),
height(height)
{
  data.resize(width*height*3);
}

void GImage::setPixel(unsigned int x, unsigned int y, std::array<unsigned int, 3> color)
{
  if(x>=width)return;
  if(y>=height)return;
  unsigned int offset = 3*(y*width+x);
  data[offset+0]=static_cast<unsigned char>(color[0]);
  data[offset+1]=static_cast<unsigned char>(color[1]);
  data[offset+2]=static_cast<unsigned char>(color[2]);
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
  std::ofstream os(filename);
  if(!os) return;
  unsigned int size = 3*width*height;
  os<<"P3"<<std::endl<<width<<" "<<height<<" 255"<<std::endl;
  for(unsigned int i=0;i<size;i++)
  {
    os<<static_cast<unsigned int>(data[i]);
    if( i%20) os<<' ';
    else os<<std::endl;
  }
  os<<std::endl;
}
