//!  A global constraint class
/*!
  Global constraint to determine minimum distance between tasks acording to a confidence percent and stochastic events.
  This is essentially a linear greater or equal constraint modified to use a distribution function instead of a common sum.
  Makes sure a level of confidence is respected to ensure robustness of a solution.

  $f$ is the confidence ratio ex: 0.95\\
  $P$ Probability function\\
  $p_i$ total duration of breakdowns.\\
  $k_i$ is the minimum distance between a task the the setup after.\\
  $n$ is the cardinality of $k$.\\
  \begin{align}
    \hspace{1cm}& \prod_{i=1}^{n}P(p_i \geq k_i) \geq f\\
    \rm{\equiv } \hspace{1cm}& \sum_{i=1}^{n}\log P(p_i \geq k_i) \geq \log f
  \end{align}


*/
#include "distribution.h"
#include "MultiPoissonDistribution.h"
#include "distributionLinearGreaterEqualConstraint.h"
#include <chuffed/core/propagator.h>
#include <iostream>
#include "chuffed_utils.h"
#include <math.h>
#include <algorithm>

#define  EPSILON 0.0000001

using namespace std;

class linearGreaterEqualDistributionConstraintPropagator : Propagator { // A_1 + A_2 + ... + A_n >= B
    vec<IntVar *> distance;

    distribution* distributionUsed;

    double confidencePercent;
    double logConfidence;

public :
    //! Constraint object constructor
    /*!
      \param distance : variable representing minimum distance
      \param confidencePercent : A confidence value to respect in a given solution. ex: 95% means that we are 95% sure the solution will stay stable in reality.
      \param pEventDuration : is an array containing all the durations of stochastic events when it happens. It is a 2D array concatenated into a 1D array. Useful in case of multiple different stochastic events are possible.
      \param pEventMeanOccurrence : the mean occurrence of any stochastic event for a given distance variable. This parameter will be used as the lambda of a distribution (Poisson in this case). Must equal the size of pEventDuration.
    */
    linearGreaterEqualDistributionConstraintPropagator(vec<IntVar *> &distance, const float confidencePercent,
                                                       distribution *distributionUsed) :
            distance(distance),
            confidencePercent(confidencePercent),
            distributionUsed(distributionUsed){
        assert(confidencePercent <= 1.0 && confidencePercent >= 0.0);

        for (size_t i = 0; i < distance.size(); ++i) {
            distance[i]->attach(this, i, EVENT_L);
        }
        logConfidence = log(confidencePercent);
    }

    //! Constraint propagator
    /*! Propagates the constraint
     * \return is the constraint satisfied
    */
    bool propagate() { //TODO:penalizer les dist du debut plus que la fin
        double sumLogMaxChance = 0;
        vec<Lit> maxLiterals;
        vec<double> logMaxChance;
        maxLiterals.reserve(distance.size());
        logMaxChance.reserve(distance.size());

        for (size_t i = 0; i < distance.size(); ++i) {
            //compute min sums
            double valMin = distributionUsed->getMin(i);
            if (valMin + EPSILON < confidencePercent) { // the percent can't be below the confidenceLevel otherwise the product will be lesser than the confidenceLevel
                vec<Lit> literals;
                int minimum = distributionUsed->calculateQuantile(i,confidencePercent);
                literals.push(distance[i]->getMinLit());//
                distance[i]->setMin(minimum+1, Reason_new(literals));
                valMin = distributionUsed->getMin(i);
            }

            //compute max sums
            double valMax = distributionUsed->getMax(i);
            double logValMax = log(valMax);
            logMaxChance.push(logValMax);
            sumLogMaxChance += logMaxChance[i];
            maxLiterals.push(distance[i]->getMaxLit());
            if(valMax < confidencePercent){

                //code to generate distributed nogoods. Use at your own risk.
//                double distribute = logConfidence - sumLogMaxChance;
//                const double toDistributePerVar = distribute / i;
//                for (size_t d = 0; d < i && distribute > 0; ++d) {
//                    double distributedPercent = min(1.0 - EPSILON,exp(logMaxChance[d] + min(distribute, toDistributePerVar)));
//                    double quantileDistributedDistance = distributionUsed->calculateQuantile(d,distributedPercent);
//
//                    maxLiterals[d] = getNegLeqLit(distance[d], distributionUsed->addMaxDistance(quantileDistributedDistance, d));
//                    distribute -= distributedPercent;
//                }
                sat.confl = Reason_new(maxLiterals);
                return false;
            }

            //if the maximum percent is lesser than the confidence percent then the schedule is risky
            if (sumLogMaxChance < logConfidence) {
                sat.confl = Reason_new(maxLiterals);
                return false;
            }
        }

        for (size_t i = 0; i < distance.size(); ++i) {
            const double logMinPossibleValueOfIndex = logConfidence - (sumLogMaxChance - logMaxChance[i]);
            double minimumPossibleValueOfIndex =  exp(logMinPossibleValueOfIndex);

            if (minimumPossibleValueOfIndex>=1){ //quantile overflows if value is 1
                minimumPossibleValueOfIndex -= EPSILON;
            }

            int minimumDistance = distributionUsed->calculateQuantile(i,minimumPossibleValueOfIndex);

            if (minimumDistance > distance[i]->getMax()) {
                vec<Lit> literals(maxLiterals);
                literals.push(getNegLeqLit(distance[i], minimumDistance - 1));//getNegLeqLit(distance[i], intMinimumDistance - 1)
                sat.confl = Reason_new(literals);
                return false;
            }

            //propagate distance[i] lower bound
            if (minimumDistance > distance[i]->getMin()) {
                vec<Lit> literals(maxLiterals);
                literals.reserve(i+1);
                //add only previous literals since others didn't contribute to propagation
                literals.push(distance[i]->getMinLit());//getNegLeqLit(distance[i], intMinimumDistance)
                distance[i]->setMin(minimumDistance, Reason_new(literals));
            }
        }
        return true;
    }
};

//! Creates the constraint
/*!
 * Creates the constraint by calling its constructor. This will create a memory leak... but chuffed works this way.
  \param distance : variable representing minimum distance
  \param confidencePercent : A confidence value to respect in a given solution. ex: 95% means that we are 95% sure the solution will stay stable in reality. Integer percent to be converted to a float [0.0,1.0].
  \param eventDuration : is an array containing all the durations of stochastic events when it happens. It is a 2D array concatenated into a 1D array. Useful in case of multiple different stochastic events are possible.
  \param eventMeanOccurrence : the mean occurrence of any stochastic event for a given distance variable. This parameter will be used as the lambda of a distribution (Poisson in this case). Must equal the size of pEventDuration. Since there is no float, this in the float value multiplied by 100 in an integer.
  \return 0
 */
int linearGreaterEqualDistributionConstraint(vec<IntVar *> &distance, const int confidencePercent,
                                             const vec<int> &eventDuration,
                                             const vec<int> &eventMeanOccurrence)
{
    distribution* poissonDistribution =  new MultiPoissonDistribution(distance,eventDuration,eventMeanOccurrence);
    new linearGreaterEqualDistributionConstraintPropagator(distance, confidencePercent / 100.0, poissonDistribution);

    return 0;
}
