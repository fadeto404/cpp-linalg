//
// Created by rns on 6/17/21.
//

#ifndef CPP_LINALG_UTIL_H
#define CPP_LINALG_UTIL_H

namespace cla
{
    template<class T>
    T max(T a, T b)
    {
        return a > b ? a : b;
    }

    template<class T>
    T min(T a, T b)
    {
        return a < b ? a : b;
    }
}


#endif //CPP_LINALG_UTIL_H
