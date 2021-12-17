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
#include <optional>
#include <cassert>


class hittable_list : public hittable  {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { objects.clear(); }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }
        pair<point3, point3> get_bounding_box() const override;
        void compile(int dim = 0) override;
        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    private:
        bool is_leaf = true;
        std::vector<shared_ptr<hittable>> objects;
        std::optional<pair<point3, point3>> bounding_box;
};


pair<point3, point3> hittable_list::get_bounding_box() const {
    assert(bounding_box.has_value());
    return bounding_box.value();
}

// Assumes each underlying object is a sphere
void hittable_list::compile(int dim) {
    assert(is_leaf == true);

    if(objects.size() <= 3){
        // do not split hierarchy
    } else {
        hittable_list left_split;
        hittable_list right_split;

        // get median, split in half
        vector<double> coords;
        for(shared_ptr<hittable> i : objects){
            double coord_min = i->get_bounding_box().first.e[dim];
            double coord_max = i->get_bounding_box().second.e[dim];

            coords.emplace_back((coord_min + coord_max) / 2.0);
        }
        double median = coords[coords.size() / 2];

        int left_added = 0;
        int right_added = 0;
        for(shared_ptr<hittable> i : objects){
            if(i->get_bounding_box().first.e[dim] <= median){
                left_split.add(i);
                left_added++;
            }
            if(i->get_bounding_box().second.e[dim] >= median){
                right_split.add(i);
                right_added++;
            }
        }

        int n = objects.size();
        // if the split is good, we make an internal node
        if(2 * (left_added + right_added) < 3 * n){
            is_leaf = false;
            left_split.compile((dim + 1) % 3);
            right_split.compile((dim + 1) % 3);

            objects = {make_shared<hittable_list>(left_split), make_shared<hittable_list>(right_split)};
        }
    }

    bool first = true;
    point3 mins;
    point3 maxs;
    for(shared_ptr<hittable> i : objects){
        auto i_box = i->get_bounding_box();
        if(first){
            mins = i_box.first;
            maxs = i_box.second;
            first = false;
        } else {
            mins.e[0] = min(mins.e[0], i_box.first.e[0]);
            mins.e[1] = min(mins.e[1], i_box.first.e[1]);
            mins.e[2] = min(mins.e[2], i_box.first.e[2]);
            maxs.e[0] = min(maxs.e[0], i_box.second.e[0]);
            maxs.e[1] = min(maxs.e[1], i_box.second.e[1]);
            maxs.e[2] = min(maxs.e[2], i_box.second.e[2]);
        }
    }

    bounding_box = {mins, maxs};
}

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const 
{
    // Do The algorithm on page 69 of https://www.cs.rochester.edu/courses/572/fall2021/decks/lect20-ray-tracing.pdf
    hit_record temp_rec;
    auto hit_anything = false;
    auto closest_so_far = t_max;

    if (!hit_the_box())
    {
        // TODO: write the code to check if the box is hitted.
        return hit_anything;
    }

    if (is_leaf)
    {
        for (const auto &object : objects)
        {
            if (object->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        return hit_anything;
    }

    // go through the two branches of the tree
    for (const auto &list : objects)
    {
        // go through all the objects
        for (const auto &object : list.objects)
        {
            if (object->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;

                rec = temp_rec;
            }
        }
    }

    return hit_anything;
}

#endif
