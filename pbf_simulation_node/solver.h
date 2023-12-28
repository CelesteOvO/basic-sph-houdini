//
// Created by LiYifan on 2023/12/27.
//

#ifndef FLUID_HOUDINI_TEST_SOLVER_H
#define FLUID_HOUDINI_TEST_SOLVER_H

#include <SIM/SIM_SingleSolver.h>
#include <SIM/SIM_OptionsUser.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_Utils.h>

class PBFSolver : public SIM_SingleSolver, public SIM_OptionsUser {
public:
private:
    explicit PBFSolver(const SIM_DataFactory *factory) : SIM_SingleSolver(factory), SIM_OptionsUser(this) {}
    ~PBFSolver() override = default;
    SIM_Result solveSingleObjectSubclass(SIM_Engine &engine, SIM_Object &object, SIM_ObjectArray &feedbacktoobjects, const SIM_Time &timestep, bool newobject) override;
    static const SIM_DopDescription *GetDescription();

    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(PBFSolver, SIM_SingleSolver, "PBF Solver Description", GetDescription());

    // ==================== Custom Field ====================
protected:
    void init(SIM_Object &obj);
    void solve(SIM_Object &obj, const SIM_Time &dt);

    struct Neighbor
    {
        Neighbor(long i, const UT_Vector3& x, double d);

        long	id;
        fpreal	dist;
        UT_Vector3	xij{};
    };
};

/// A simple Logger
#include <iostream>
#include <vector>
#include <string>

static struct Logger
{
    SIM_Solver::SIM_Result report()
    {
        for (const auto &l: log) std::cout << l << '\n';
        return log.empty() ? SIM_Solver::SIM_SOLVER_SUCCESS : SIM_Solver::SIM_SOLVER_FAIL;
    }
    void reset() { log.clear(); }
    void error_nullptr(const std::string &type) { log.emplace_back("NullPointerError::" + type); }

    std::vector<std::string> log;
} Log;

#endif //FLUID_HOUDINI_TEST_SOLVER_H
