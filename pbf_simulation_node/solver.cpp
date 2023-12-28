//
// Created by LiYifan on 2023/12/27.
//

#include "solver.h"
#include "sph_kernel.h"
#include "SpatialGrid.h"

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

const SIM_DopDescription *PBFSolver::GetDescription()
{
    static std::array<PRM_Template, 2> PRMS{
            PRM_Template()
    };
    static SIM_DopDescription DESC(
            true,
            "pbf_solverCeleste",
            "PBF SolverCeleste",
            "PBFSolverCeleste",
            classname(),
            PRMS.data()
    );
    return &DESC;
}

SIM_Solver::SIM_Result
PBFSolver::solveSingleObjectSubclass(SIM_Engine &engine, SIM_Object &object, SIM_ObjectArray &feedbacktoobjects,
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

void PBFSolver::init(SIM_Object &obj) {
    /*std::cout << "Initializing simulator..." << std::endl;*/

    SIM_GeometryCopy *geo;
    geo = SIM_DATA_CREATE(obj, "Geometry", SIM_GeometryCopy,
                          SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE);
    if (!geo)
        Log.error_nullptr("INIT::SIM_GeometryCopy");

    {
        SIM_GeometryAutoWriteLock lock(geo);
        GU_Detail &gdp = lock.getGdp();

        GA_RWAttributeRef mass_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS), 1, GA_Defaults(0));
        GA_RWAttributeRef vel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, gdp.getStdAttributeName(GEO_ATTRIBUTE_VELOCITY),
                                                      3, GA_Defaults(0));

        GA_RWAttributeRef predict_pos_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "predict_pos", 3, GA_Defaults(0));

        GA_RWAttributeRef lambda_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "lambda", 1, GA_Defaults(0));
        GA_RWAttributeRef delta_p_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "delta_p", 3, GA_Defaults(0));

        GA_RWAttributeRef density_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "density", 1, GA_Defaults(0));
        GA_RWAttributeRef pressure_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pressure", 1, GA_Defaults(0));

        GA_RWAttributeRef external_accel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "eA", 3, GA_Defaults(0));
        GA_RWAttributeRef pressure_accel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pA", 3, GA_Defaults(0));

        GA_RWAttributeRef neighbor_num_ref = gdp.addIntTuple(GA_ATTRIB_POINT, "nNum", 1, GA_Defaults(0));


        GA_RWHandleF masshandle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));
        GA_RWHandleV3 pos_handle = gdp.findPointAttribute("P");

        fpreal target_spacing = 0.029;
        fpreal relative_kernel_radius = 1.7;
        fpreal kernel_radius = target_spacing * relative_kernel_radius;

        fpreal target_density = 1000.0;

        fpreal h = kernel_radius;
        fpreal radius = 0.029;
        UT_Vector3 volumeMin = UT_Vector3 (-0.5, 0, -0.5);
        UT_Vector3 volumeMax = UT_Vector3 (0.5, 1.5, 0.5);

        long long particle_size = gdp.getNumPoints();

        /// 1.build Neighbors
        std::vector<std::vector<Neighbor> >	_temp_neighbors;
        // Reserve space and initialize neighbors' data
        _temp_neighbors.clear();
        _temp_neighbors.resize(gdp.getNumPoints());
        // Init spatial grid
        UT_Vector3 borders(h*2.0, h*2.0, h*2.0);
        UT_Vector3 gridMin = volumeMin;
        gridMin -= borders;
        UT_Vector3 gridMax = volumeMax;
        gridMax += borders;
        SpatialGrid<long> grid(h, gridMin, gridMax);
        // Insert particles into grid
        for (long p = 0; p < particle_size; ++p)
        {
            GA_Offset offset = gdp.pointOffset(p);
            UT_Vector3 pos = pos_handle.get(offset);
            grid.insert(p, pos);
        }
        // Use grid to retrieve neighbor particles
        double h2 = h*h;
        std::vector<long*> nearbyParticles;

        for (long p=0; p<particle_size; ++p)
        {
            int neighbor_num = 0;
            GA_Offset offset = gdp.pointOffset(p);
            UT_Vector3 pos = pos_handle.get(offset);
            // Get nearby particles
            grid.getElements(pos, h, nearbyParticles);

            // Find particles that are within smoothing radius
            _temp_neighbors[p].reserve(50);
            for (auto & nearbyParticle : nearbyParticles)
            {
                long nID = *nearbyParticle;
                GA_Offset nearbyOffset = gdp.pointOffset(nID);
                UT_Vector3 nearbyPos = pos_handle.get(nearbyOffset);

                // Skip current particle
                if (nID==p)
                    continue;

                UT_Vector3 xij = pos - nearbyPos;

                // Check if distance is lower than smoothing radius
                double dist2 = xij.dot(xij);
                if (dist2 < h2)
                {
                    // Yup! Add the particle to the neighbors list along with
                    // some precomputed information
                    _temp_neighbors[p].emplace_back(nID, xij, sqrt(dist2));
                    neighbor_num++;
                }
            }
        }

        StdKernel poly6(kernel_radius);
        fpreal max_number_density = 0;
        for(long long i = 0; i < particle_size; i++) {
            GA_Offset offset = gdp.pointOffset(i);
            fpreal sum = poly6(0.0);
            for (const auto & neighbor : _temp_neighbors[i]) {
                // Add contribution
                sum += poly6(neighbor.dist);
            }
            max_number_density = std::max(max_number_density, sum);
        }

        if(max_number_density > 0)
        {

            fpreal mass = target_density / max_number_density;
            for(long long i = 0; i < particle_size; i++) {
                GA_Offset offset = gdp.pointOffset(i);
                masshandle.set(offset, mass);
            }
        }

        /// 2.compute density
        GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));
        GA_RWHandleF density_handle = gdp.findPointAttribute("density");
        GA_RWHandleF pressure_handle = gdp.findPointAttribute("pressure");
        GA_RWHandleV3 external_accel_handle = gdp.findPointAttribute("eA");
        GA_RWHandleV3 predict_pos_handle = gdp.findPointAttribute("predict_pos");

        for(long long i = 0; i < particle_size; i++) {
            GA_Offset offset = gdp.pointOffset(i);
            predict_pos_handle.set(i, pos_handle.get(offset));
        }

        for (long p=0; p<particle_size; ++p)
        {
            GA_Offset offset = gdp.pointOffset(p);

            // Reinitialize particle properties
            density_handle.set(offset, 0.0);
            external_accel_handle.set(offset, UT_Vector3(0,0,0));
            pressure_handle.set(offset, 0.0);

            fpreal particle_density = mass_handle.get(offset) * poly6(0.0);
            for (const auto & neighbor : _temp_neighbors[p])
            {
                // Add contribution
                particle_density += mass_handle.get(neighbor.id) * poly6(neighbor.dist);
            }
            density_handle.set(offset, particle_density);
        }
    }
}

void PBFSolver::solve(SIM_Object &obj, const SIM_Time &timestep) {
    /*std::cout << "Solving..." << std::endl;*/
    SIM_GeometryCopy *geo;
    geo = SIM_DATA_CREATE(obj, "Geometry", SIM_GeometryCopy,
                          SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE);
    if (!geo)
        Log.error_nullptr("SOLVE::SIM_GeometryCopy");

    {
        SIM_GeometryAutoWriteLock lock(geo);
        GU_Detail &gdp = lock.getGdp();

        GA_RWHandleV3 pos_handle = gdp.findPointAttribute("P");
        GA_RWHandleV3 vel_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_VELOCITY));
        GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));
        GA_RWHandleV3 predict_pos_handle = gdp.findPointAttribute("predict_pos");
        GA_RWHandleF lambda_handle = gdp.findPointAttribute("lambda");
        GA_RWHandleV3 delta_p_handle = gdp.findPointAttribute("delta_p");
        GA_RWHandleF density_handle = gdp.findPointAttribute("density");
        GA_RWHandleF pressure_handle = gdp.findPointAttribute("pressure");
        GA_RWHandleV3 external_accel_handle = gdp.findPointAttribute("eA");
        GA_RWHandleV3 pressure_accel_handle = gdp.findPointAttribute("pA");

        fpreal target_spacing = 0.029;
        fpreal relative_kernel_radius = 1.7;
        fpreal kernel_radius = target_spacing * relative_kernel_radius;

        fpreal target_density = 1000.0;

        fpreal _h = kernel_radius;
        fpreal radius = 0.029;
        UT_Vector3 _volumeMin = UT_Vector3 (-0.5, 0, -0.5);
        UT_Vector3 _volumeMax = UT_Vector3 (0.5, 1.5, 0.5);
        fpreal dt = 0.005;


        long long particle_size = gdp.getNumPoints();

        /// 1. update position and velocity
        for(long long i = 0; i < particle_size; i++) {
            GA_Offset offset = gdp.pointOffset(i);

            UT_Vector3 velocity = (predict_pos_handle.get(offset) - pos_handle.get(offset)) / dt;
            UT_Vector3 position = predict_pos_handle.get(offset);

            if(position.y() < 0) {
                position.y() = 0;
                velocity.y() = 0;
            }

            vel_handle.set(offset, velocity);
            pos_handle.set(offset, position);
        }

        /*for(long long i = 0; i < particle_size; i++) {
            GA_Offset offset = gdp.pointOffset(i);

            UT_Vector3 position = predict_pos_handle.get(offset);
            pos_handle.set(offset, position);
        }*/

        /// 2. apply force and predict position
        for(long long i = 0; i < particle_size; i++) {
            GA_Offset offset = gdp.pointOffset(i);
            predict_pos_handle.set(offset, pos_handle.get(offset));

            UT_Vector3 external_accel = UT_Vector3(0,-9.81,0);
            external_accel_handle.set(offset, external_accel);

            UT_Vector3 velocity = vel_handle.get(offset);
            velocity += external_accel * dt;
            vel_handle.set(offset, velocity);

            UT_Vector3 position = predict_pos_handle.get(offset);
            position += velocity * dt;
            predict_pos_handle.set(offset, position);
        }

        /// 3.solve density constraints
        int constraint_solver_iterations = 5;
        for(int m = 0; m < constraint_solver_iterations; ++m)
        {
            /// 3.1.build Neighbors
            std::vector<std::vector<Neighbor> >	_neighbors;
            // Reserve space and initialize neighbors' data
            _neighbors.clear();
            _neighbors.resize(particle_size);
            // Init spatial grid
            UT_Vector3 borders(_h*2.0, _h*2.0, _h*2.0);
            UT_Vector3 gridMin = _volumeMin;
            gridMin -= borders;
            UT_Vector3 gridMax = _volumeMax;
            gridMax += borders;
            SpatialGrid<long> grid(_h, gridMin, gridMax);
            // Insert particles into grid
            for (long p = 0; p < particle_size; ++p)
            {
                GA_Offset offset = gdp.pointOffset(p);
                UT_Vector3 pos = predict_pos_handle.get(offset);
                grid.insert(p, pos);
            }
            // Use grid to retrieve neighbor particles
            double h2 = _h*_h;
            std::vector<long*> nearbyParticles;
            for (long p=0; p<particle_size; ++p)
            {
                GA_Offset offset = gdp.pointOffset(p);
                UT_Vector3 pos = predict_pos_handle.get(offset);
                // Get nearby particles
                grid.getElements(pos, _h, nearbyParticles);

                // Find particles that are within smoothing radius
                _neighbors[p].reserve(50);
                for (auto & nearbyParticle : nearbyParticles)
                {
                    long nID = *nearbyParticle;
                    GA_Offset nearbyOffset = gdp.pointOffset(nID);
                    UT_Vector3 nearbyPos = predict_pos_handle.get(nearbyOffset);

                    // Skip current particle
                    if (nID==p)
                        continue;

                    UT_Vector3 xij = pos - nearbyPos;

                    // Check if distance is lower than smoothing radius
                    double dist2 = xij.dot(xij);
                    if (dist2 < h2)
                    {
                        // Yup! Add the particle to the neighbors list along with
                        // some precomputed informations
                        _neighbors[p].emplace_back(nID, xij, sqrt(dist2));
                    }
                }
            }

            /// 3.2.compute density
            StdKernel poly6(kernel_radius);
            for (long p=0; p<particle_size; ++p)
            {
                GA_Offset offset = gdp.pointOffset(p);

                fpreal particle_density = mass_handle.get(offset) * poly6(0.0);
                for (const auto & neighbor : _neighbors[p])
                {
                    // Add contribution
                    fpreal dist = (predict_pos_handle.get(offset) - predict_pos_handle.get(gdp.pointOffset(neighbor.id))).length();
                    particle_density += mass_handle.get(neighbor.id) * poly6(dist);
                }
                density_handle.set(offset, particle_density);
            }

            /// 3.3.compute lambda
            for(long long i = 0; i < particle_size; i++) {
                GA_Offset offset = gdp.pointOffset(i);

                fpreal d_i = density_handle.get(offset);
                fpreal C_i = d_i / target_density - 1.0;

                if(C_i > 0) {
                    fpreal sum_grad_C_i_squared = 0;
                    UT_Vector3 grad_C_i = UT_Vector3(0,0,0);

                    for (const auto & neighbor : _neighbors[i])
                    {
                        GA_Offset neighbor_offset = gdp.pointOffset(neighbor.id);
                        const auto p_i = predict_pos_handle.get(offset);
                        const auto p_j = predict_pos_handle.get(neighbor_offset);
                        const UT_Vector3 grad_C_j = -(mass_handle.get(offset) / target_density) * poly6.gradient(p_i - p_j);

                        sum_grad_C_i_squared += grad_C_j.dot(grad_C_j);
                        grad_C_i -= grad_C_j;
                    }

                    sum_grad_C_i_squared += grad_C_i.dot(grad_C_i);

                    fpreal lambda = -C_i / (sum_grad_C_i_squared + 1e-6);
                    lambda_handle.set(offset, lambda);
                }else{
                    lambda_handle.set(offset, 0);
                }
            }

            for(long long i = 0; i < particle_size; i++) {
                GA_Offset offset = gdp.pointOffset(i);

                fpreal lambda_i = lambda_handle.get(offset);

                auto k_corr = mass_handle.get(offset) * 1.0e-04;
                auto n_corr = 4.0;
                auto q_corr = 0.1;

                UT_Vector3 delta_p_i = UT_Vector3(0,0,0);

                for (const auto & neighbor : _neighbors[i])
                {
                    GA_Offset neighbor_offset = gdp.pointOffset(neighbor.id);
                    const auto lambda_j = lambda_handle.get(neighbor_offset);
                    const auto p_i = predict_pos_handle.get(offset);
                    const auto p_j = predict_pos_handle.get(neighbor_offset);

                    const auto w_corr = poly6(q_corr * radius);
                    const auto ratio = poly6((p_i - p_j).length()) / w_corr;
                    const auto s_corr = -k_corr * pow(ratio, n_corr);

                    const UT_Vector3 grad_C_j = -(mass_handle.get(offset) / target_density) * poly6.gradient(p_i - p_j);

                    delta_p_i -= (lambda_i + lambda_j + s_corr) * grad_C_j;
                }
                delta_p_handle.set(offset, delta_p_i);
            }

            /// 3.4.update position
            for(long long i = 0; i < particle_size; i++) {
                GA_Offset offset = gdp.pointOffset(i);

                UT_Vector3 p_to_write = predict_pos_handle.get(offset);
                p_to_write -= delta_p_handle.get(offset);
                predict_pos_handle.set(offset, p_to_write);
            }
        }
    }
}

PBFSolver::Neighbor::Neighbor(long i, const UT_Vector3& x, double d) {
    id = i;
    xij = x;
    dist = d;
}


