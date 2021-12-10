#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H
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

#include "rtweekend.h"

#include "hittable.h"

#include <memory>
#include <vector>


class hittable_list : public hittable  {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }
        pair<vec3, vec3> get_bounding_box() const override;
        void compile() const override;
        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    private:
        bool is_leaf = true;
        std::vector<shared_ptr<hittable>> objects;
};


pair<vec3, vec3> hittable_list::get_bounding_box() const {
    // TODO
    // For each object, get the bounding boxes
    // Take the minimum x and y coordinate and maximum x and y coordinate as bounding box
}

// Assumes each underlying object is a sphere
void hittable_list::compile() const {
    // TODO
    assert(is_leaf == true)
    // If there is <= 3 objects
    //     make it a leaf.
    // Otherwise, make a split.
    //     Sort elements by x coordinate, get a_1, a_2, ... a_n
    //     Put a_1, a_2, ...., a_n/2 in one half
    //     Put a_n/2+1, ... a_n in other half
    //     set is_leaf = true
}

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // TODO
    // If is_leaf
        // hit_record temp_rec;
        // auto hit_anything = false;
        // auto closest_so_far = t_max;

        // for (const auto& object : objects) {
        //     if (object->hit(r, t_min, closest_so_far, temp_rec)) {
        //         hit_anything = true;
        //         closest_so_far = temp_rec.t;
        //         rec = temp_rec;
        //     }
        // }
        //
        // return hit_anything;
    // Else
        // Do The algorithm on page 69 of https://www.cs.rochester.edu/courses/572/fall2021/decks/lect20-ray-tracing.pdf
}


#endif
