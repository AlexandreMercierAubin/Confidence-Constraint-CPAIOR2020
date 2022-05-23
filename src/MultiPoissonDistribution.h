//
// Created by alexandre mercier-aubin on 4/3/20.
//

#ifndef STOCHASTIC_CUMULATIVE_MULTIPOISSONDISTRIBUTION_H
#define STOCHASTIC_CUMULATIVE_MULTIPOISSONDISTRIBUTION_H
#include <chuffed/vars/int-var.h>
#include <vector>
#include "distribution.h"
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/poisson.hpp>
using namespace std;
using namespace boost::math;
class MultiPoissonDistribution: public distribution {
private:
    vec<IntVar*> variables;
    vector<int> eventDuration;
    vector<float> eventMeanOccurrence;
    int eventTypesCount;
    vector<size_t> distributionMap;
    vector<poisson_distribution<>> poissonDistributions;

    double getValue(const double k, const size_t index) const;
    double calculateLambda(const size_t index) const;
    double ponderateK(const double k, const size_t index) const;
public:
    MultiPoissonDistribution(vec<IntVar*> &variables, const vec<int> &pEventDuration,
                             const vec<int> &pEventMeanOccurrence);
    double getMax(const size_t index) const override;
    double getMin(const size_t index) const override;
    double calculateQuantile(const size_t index, const double percent) const override;
    double addMaxDistance(double distanceToAdd, const size_t index) const override;
};


#endif //STOCHASTIC_CUMULATIVE_MULTIPOISSONDISTRIBUTION_H
