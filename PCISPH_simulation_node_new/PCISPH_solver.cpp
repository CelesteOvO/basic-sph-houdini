//
// Created by LiYifan on 2023/12/25.
//

#include "PCISPH_solver.h"
#include "include/BoundingBox.h"
#include "include/PciSphSolver.h"

#include <PRM/PRM_Name.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Shared.h>
#include <PRM/PRM_Default.h>

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_ObjectArray.h>
#include <SIM/SIM_GeometryCopy.h>

#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_Primitive.h>

#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_ParallelUtil.h>
#include <UT/UT_ParallelPipeline.h>

#include <array>
#include <iostream>

const SIM_DopDescription *PCISPHSolver::GetDescription()
{
    static std::array<PRM_Template, 1> PRMS{
            PRM_Template()
    };

    static SIM_DopDescription DESC(true,
                                   "pcisph_solver",
                                   "PCISPH Solver",
                                   "PCISPHSolver",
                                   classname(),
                                   PRMS.data());
    return &DESC;
}

SIM_Solver::SIM_Result
PCISPHSolver::solveSingleObjectSubclass(SIM_Engine &engine, SIM_Object &object, SIM_ObjectArray &feedbacktoobjects,
                                        const SIM_Time &timestep, bool newobject) {
    Log.reset();
    static bool NeedReBuild = true;
    if (NeedReBuild || newobject) {
        init(object);
        NeedReBuild = false;
    }else{
        solve(object, timestep);
    }
    return Log.report();
}

void PCISPHSolver::init(SIM_Object &obj) const {
    SIM_GeometryCopy *geo;
    geo = SIM_DATA_CREATE(obj, "Geometry", SIM_GeometryCopy,
                          SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE);
    if (!geo)
        Log.error_nullptr("INIT::SIM_GeometryCopy");
}

void PCISPHSolver::solve(SIM_Object &obj, const SIM_Time &dt) const
{
    double targetSpacing = 0.02;
    double targetDensity = 1000;
    double relativeKernelRadius = 1.8;

    UT_Vector3 lowerCorner, upperCorner;
    BoundingBox domain(lowerCorner, upperCorner);

    // Build solver
    PciSphSolver solver = PciSphSolver(targetDensity, targetSpacing, relativeKernelRadius);

    // Build collider

}