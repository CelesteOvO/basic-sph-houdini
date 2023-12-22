//
// Created by LiYifan on 2023/12/4.
//
#include <UT/UT_DSOVersion.h> // Very Important!!! Include this!!!

#include "solver.h"

void initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(PCISPHSolver);
}
