#ifndef D6_SHAPE_RENDERER_HH
#define D6_SHAPE_RENDERER_HH

#include "utils/alg.hh"

namespace d6 {

class LevelSet ;

class ShapeRenderer
{

public:
    void init() ;

    void draw(const LevelSet &ls, const Vec &box ) const ;

private:

};

} //d6

#endif
