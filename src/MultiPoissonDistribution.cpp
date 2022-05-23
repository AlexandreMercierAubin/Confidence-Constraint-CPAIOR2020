//
// Created by alexandre mercier-aubin on 4/3/20.
//

#include "MultiPoissonDistribution.h"
#include <boost/math/distributions/poisson.hpp>

using namespace boost::math;

//! Constructor
/*! Inits variables
*/
MultiPoissonDistribution::MultiPoissonDistribution(vec<IntVar *> &variables, const vec<int> &pEventDuration,
                                                   const vec<int> &pEventMeanOccurrence) : variables(variables) {
    for (size_t i = 0; i < pEventDuration.size(); ++i) {
        eventDuration.push_back(pEventDuration[i]);
        eventMeanOccurrence.push_back(pEventMeanOccurrence[i] / 100.0);
    }
    eventTypesCount = pEventMeanOccurrence.size() / variables.size();
    for(size_t i = 0; i<variables.size(); ++i){

        double lambda = calculateLambda(i);
        if (lambda > 0){
            distributionMap.push_back(poissonDistributions.size());
            poisson_distribution<> p(lambda);
            poissonDistributions.push_back(p);
        }
        else{
            distributionMap.push_back(-1);
        }
    }


}

//! Gets maximum value of the domain of the variable at index
/*!
 * \return percent
*/
double MultiPoissonDistribution::getMax(const size_t index) const {
    return getValue(variables[index]->getMax(), index);
}

//! Gets minimum value of the domain of the variable at index
/*!
 * \return percent
*/
double MultiPoissonDistribution::getMin(const size_t index) const {
    return getValue(variables[index]->getMin(), index);
}

//! Get a value of variable at index, k is a distance
/*! Calculates the distance according to a percent and event particularities.
 * This is the propability function of this distribution.
 * \return percent
*/
double MultiPoissonDistribution::getValue(const double k, const size_t index) const {
    //This is ok since we use floats for percents. Otherwise, we could add the probability of exact equality
    double result = 1.0;
    size_t pos = distributionMap[index];
    if(pos != -1) {
        result = boost::math::cdf(poissonDistributions[pos], k);
    }
    return result;
}

//! Quantile function
/*! Calculates the distance according to a percent and event particularities
 * \return distance
*/
double MultiPoissonDistribution::calculateQuantile(size_t index, double percent) const {
    size_t pos = distributionMap[index];
//    if (percent >= 1.0){
//        return variables[index]->getMax();
//    }
//    if (percent < 0.0){
//        return variables[index]->getMin();
//    }
    double distance = 0.0;
    if (pos !=-1) {
        distance = quantile(poissonDistributions[pos], percent);
    }
    return distance;
}

//! Calculates the ponderated lambda
/*!
 * \return ponderated lambda parameter
*/
double MultiPoissonDistribution::calculateLambda(const size_t index) const {
    double lambda = 0.0;
    for (size_t i = 0; i < eventTypesCount; ++i) {
        size_t pos = index + i * variables.size();
        lambda += eventMeanOccurrence[pos]*eventDuration[pos];
    }
    return lambda;
}

//! Adds a value to the upper bound of variable i
/*!
 * \return the sum
*/
double MultiPoissonDistribution::addMaxDistance(double distanceToAdd, const size_t index) const {
    return variables[index]->getMax() + distanceToAdd;
}
