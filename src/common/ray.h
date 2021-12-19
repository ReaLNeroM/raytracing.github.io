#ifndef RAY_H
#define RAY_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "vec3.h"


class ray {
    public:
        ray() {}
        ray(const point3& origin, const vec3& direction)
            : orig(origin), dir(direction), tm(0)
        {
            invdir = vec3(1. / dir.e[0], 1. / dir.e[1], 1. / dir.e[2]);
            sign[0] = (invdir.x() < 0);
            sign[1] = (invdir.y() < 0);
            sign[2] = (invdir.z() < 0);
        }

        ray(const point3& origin, const vec3& direction, double time)
            : orig(origin), dir(direction), tm(time)
        {
            invdir = vec3(1. / dir.e[0], 1. / dir.e[1], 1. / dir.e[2]);
            sign[0] = (invdir.x() < 0);
            sign[1] = (invdir.y() < 0);
            sign[2] = (invdir.z() < 0);
        }

        point3 origin() const  { return orig; }
        vec3 direction() const { return dir; }
        double time() const    { return tm; }

        point3 at(double t) const {
            return orig + t*dir;
        }

    public:
        point3 orig;
        vec3 dir;
        vec3 invdir;
        int sign[3];
        double tm;
};

#endif
