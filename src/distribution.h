//
// Created by alexandre mercier-aubin on 4/3/20.
//

#ifndef STOCHASTIC_CUMULATIVE_DISTRIBUTION_H
#define STOCHASTIC_CUMULATIVE_DISTRIBUTION_H
#include <cstddef>
using namespace std;

class distribution{
public:
    distribution(){};
    virtual double getMax(const size_t index) const = 0;
    virtual double getMin(const size_t index) const = 0;
    virtual double calculateQuantile(const size_t index, const double percent) const = 0;
    virtual double addMaxDistance(double distanceToAdd, const size_t index) const = 0;
};

#endif //STOCHASTIC_CUMULATIVE_DISTRIBUTION_H
