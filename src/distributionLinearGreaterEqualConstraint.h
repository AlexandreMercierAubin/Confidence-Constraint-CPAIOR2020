#ifndef DISTRIBUTION_LINEAR_GREATER_EQUAL
#define DISTRIBUTION_LINEAR_GREATER_EQUAL
#include <chuffed/vars/int-var.h>

int linearGreaterEqualDistributionConstraint(vec<IntVar*> &distance, const int confidencePercent, const vec<int> &eventDuration, const vec<int> &eventMeanOccurrence);

#endif