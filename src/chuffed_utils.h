//
// Created by user on 2/18/20.
//

#ifndef STOCHASTIC_CUMULATIVE_CHUFFED_UTILS_H
#define STOCHASTIC_CUMULATIVE_CHUFFED_UTILS_H

#include <chuffed/vars/int-var.h>

Lit getNegLeqLit(IntVar *v, int val);
Lit getNegGeqLit(IntVar *v, int val);

#endif //STOCHASTIC_CUMULATIVE_CHUFFED_UTILS_H
