//
// Created by user on 2/18/20.
//
#include "chuffed_utils.h"

Lit getNegLeqLit(IntVar *v, int val){
    return (INT_VAR_LL == v->getType() ? v->getMaxLit() : v->getLit(val + 1, 2));
}
Lit getNegGeqLit(IntVar *v, int val){
    return (INT_VAR_LL == v->getType() ? v->getMaxLit() : v->getLit(val - 1, 3));
}