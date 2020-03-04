/* inclusion guard */
#ifndef __GImage_H__
#define __GImage_H__

#include <string>
#include <fstream>
#include <vector>
#include <array>

class GImage
{
public:
  GImage(unsigned int width, unsigned int height);
  void setPixel(unsigned int x,unsigned int y,std::array<unsigned int,3> color);
  std::array<unsigned int,3> getPixel(unsigned int x,unsigned int y) const;
  unsigned int getWidth() const;
  unsigned int getHeight() const;
  void save(std::string filename) const;
private:
  unsigned int width;
  unsigned int height;
  std::vector<unsigned char> data;
};


#endif /* __GImage_H__ */
