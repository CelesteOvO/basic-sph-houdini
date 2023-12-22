//
// Created by LiYifan on 2023/12/21.
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

PCISPHSolver::Neighbor::Neighbor(long i, const UT_Vector3& x, double d) {
    id = i;
    xij = x;
    dist = d;
}

void PCISPHSolver::init(SIM_Object &obj) const {
    SIM_GeometryCopy *geo;
    geo = SIM_DATA_CREATE(obj, "Geometry", SIM_GeometryCopy,
                          SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE);
    if (!geo)
        Log.error_nullptr("INIT::SIM_GeometryCopy");

    {
        SIM_GeometryAutoWriteLock lock(geo);
        GU_Detail &gdp = lock.getGdp();

        GA_RWAttributeRef vel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, gdp.getStdAttributeName(GEO_ATTRIBUTE_VELOCITY), 3, GA_Defaults(0));
        GA_RWAttributeRef mass_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "mass", 1, GA_Defaults(0));
        GA_RWAttributeRef density_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "density", 1, GA_Defaults(0));
        GA_RWAttributeRef pressure_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pressure", 1, GA_Defaults(0));

        GA_RWAttributeRef predict_position_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pP", 3, GA_Defaults(0));
        GA_RWAttributeRef predict_velocity_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pV", 3, GA_Defaults(0));
        GA_RWAttributeRef predict_density_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pDensity", 1, GA_Defaults(0));

        GA_RWAttributeRef density_error_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "Derr", 1, GA_Defaults(0));

        GA_RWAttributeRef external_force_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "eF", 3, GA_Defaults(0));
        GA_RWAttributeRef pressure_force_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pF", 3, GA_Defaults(0));

        GA_RWAttributeRef neighbor_num_ref = gdp.addIntTuple(GA_ATTRIB_POINT, "nNum", 1, GA_Defaults(0));

        //////////////////////////////////////////////////////////////////////
        fpreal radius 			= 0.029;
        fpreal target_density 	= 1000; // water density
        fpreal target_spacing 	= radius;
        fpreal relative_kernel_radius = 1.7;
        fpreal kernel_radius 		= target_spacing * relative_kernel_radius;

        UT_Vector3 volumeMin = UT_Vector3 (-3, 0, -1);
        UT_Vector3 volumeMax = UT_Vector3 (3, 3, 1);
        fpreal h = kernel_radius;
        /////////////////////////////////////////////////////////////////////

        GA_RWHandleV3 pos_handle = gdp.findPointAttribute("P");
        GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));

        long long particle_size = gdp.getNumPoints();

        std::vector<std::vector<Neighbor> >	temp_neighbor_list;
        temp_neighbor_list.resize(particle_size);
        UT_Vector3 borders(h * 2.0, h * 2.0, h * 2.0);
        UT_Vector3 gridMin = volumeMin;
        gridMin -= borders;
        UT_Vector3 gridMax = volumeMax;
        gridMax += borders;

        SpatialGrid<long> grid(h, gridMin, gridMax);
        for (long p = 0; p < particle_size; ++p)
        {
            GA_Offset offset = gdp.pointOffset(p);
            UT_Vector3 pos = pos_handle.get(offset);
            grid.insert(p, pos);
        }

        // Use grid to retrieve neighbor particles
        double h2 = h*h;
        std::vector<long*> nearbyParticles;
        for (long p = 0; p < particle_size; ++p)
        {
            GA_Offset offset = gdp.pointOffset(p);
            UT_Vector3 pos = pos_handle.get(offset);
            // Get nearby particles
            grid.getElements(pos, h, nearbyParticles);

            // Find particles that are within smoothing radius
            temp_neighbor_list[p].reserve(50);
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
                    temp_neighbor_list[p].emplace_back(nID, xij, sqrt(dist2));
                }
            }
        }

        StdKernel poly6(kernel_radius);
        fpreal max_number_density = 0;
        for (int i = 0; i < particle_size; ++i)
        {
            GA_Offset offset = gdp.pointOffset(i);

            fpreal sum = poly6(0); // self density
            for (const auto &neighbor_point_id: temp_neighbor_list[i])
            {
                sum += poly6(neighbor_point_id.dist);
            }
            max_number_density = std::max(max_number_density, sum);
        }


        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            if (max_number_density > 0)
            {
                fpreal mass = std::max((target_density / max_number_density), 0.0);
                mass_handle.set(offset, mass);
            }
            else
                throw std::runtime_error("max_number_density is zero");
        }
    }
}

void PCISPHSolver::solve(SIM_Object &obj, const SIM_Time &dt) const {
    SIM_GeometryCopy *geo;
    geo = SIM_DATA_GET(obj, "Geometry", SIM_GeometryCopy);
    if (!geo)
    {
        Log.error_nullptr("SOLVE::SIM_GeometryCopy");
        return;
    }

    {
        SIM_GeometryAutoWriteLock lock(geo);
        GU_Detail &gdp = lock.getGdp();

        GA_RWHandleV3 pos_handle = gdp.findPointAttribute("P");
        GA_RWHandleV3 vel_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_VELOCITY));
        GA_RWHandleF mass_handle = gdp.findPointAttribute("mass");
        GA_RWHandleF density_handle = gdp.findPointAttribute("density");
        GA_RWHandleF pressure_handle = gdp.findPointAttribute("pressure");

        GA_RWHandleV3 predict_position_handle = gdp.findPointAttribute("pP");
        GA_RWHandleV3 predict_velocity_handle = gdp.findPointAttribute("pV");
        GA_RWHandleF predict_density_handle = gdp.findPointAttribute("pDensity");
        GA_RWHandleF density_error_handle = gdp.findPointAttribute("Derr");

        GA_RWHandleV3 external_force_handle = gdp.findPointAttribute("eF");
        GA_RWHandleV3 pressure_force_handle = gdp.findPointAttribute("pF");

        GA_RWHandleI neighbor_num_handle = gdp.findPointAttribute("nNum");

        //////////////////////////////////////////////////////////////////////
        fpreal radius 			= 0.029;
        fpreal target_density 	= 1000; // water density
        fpreal target_spacing 	= radius;
        fpreal relative_kernel_radius = 1.7;
        fpreal kernel_radius 		= target_spacing * relative_kernel_radius;

        UT_Vector3 volumeMin = UT_Vector3 (-3, 0, -1);
        UT_Vector3 volumeMax = UT_Vector3 (3, 3, 1);
        fpreal h = kernel_radius;

        int iteration = 3;
        UT_Vector3 gravity(0, -9.81, 0);
        /////////////////////////////////////////////////////////////////////

        long long particle_size = gdp.getNumPoints();
        /// 1.build Neighbors
        std::vector<std::vector<Neighbor> >	_neighbors;
        // Reserve space and initialize neighbors' data
        _neighbors.clear();
        _neighbors.resize(particle_size);
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

        int neighbor_num = 0;

        for (long p=0; p<particle_size; ++p)
        {
            GA_Offset offset = gdp.pointOffset(p);
            UT_Vector3 pos = pos_handle.get(offset);
            // Get nearby particles
            grid.getElements(pos, h, nearbyParticles);

            // Find particles that are within smoothing radius
            _neighbors[p].reserve(50);
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
                    _neighbors[p].emplace_back(nID, xij, sqrt(dist2));
                    neighbor_num++;
                    neighbor_num_handle.set(offset, neighbor_num);
                }
            }
        }

        /// 2.compute density
        StdKernel poly6(kernel_radius);
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            fpreal sum = mass_handle(offset) * poly6(0); // self density
            for (const auto &neighbor_point_id: _neighbors[i])
            {
                sum += mass_handle(offset) * poly6(neighbor_point_id.dist);
            }
            density_handle.set(offset, sum);
        }

        /// 3.add external force
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            UT_Vector3 external_force = gravity * mass_handle.get(offset);
            external_force_handle.set(offset, external_force);
        }

        /// 4.initialize pressure and pressure force
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            pressure_handle.set(offset, 0);
            pressure_force_handle.set(offset, UT_Vector3(0, 0, 0));
            density_error_handle.set(offset, 0);
        }

        /// 5.prediction and correction
        int loop = 0;
        while(loop < iteration)
        {
            /// 5.1 predict velocity and position
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                UT_Vector3 velocity = vel_handle.get(offset);
                UT_Vector3 position = pos_handle.get(offset);

                UT_Vector3 predict_velocity = velocity + dt * external_force_handle.get(offset) / mass_handle.get(offset);
                UT_Vector3 predict_position = position + dt * predict_velocity;

                predict_velocity_handle.set(offset, predict_velocity);
                predict_position_handle.set(offset, predict_position);
            }

            /// 5.2 predict density
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                fpreal sum = mass_handle.get(offset) * poly6(0); // self density
                for (const auto &neighbor_point_id: _neighbors[i])
                {
                    fpreal predict_dist = (predict_position_handle.get(offset) - predict_position_handle.get(neighbor_point_id.id)).length();
                    sum += mass_handle.get(offset) * poly6(predict_dist);
                }
                predict_density_handle.set(offset, sum);
            }

            /// 5.3 update pressure
            fpreal delta = _compute_delta(gdp);
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                fpreal predict_density = predict_density_handle.get(offset);
                fpreal density_error = predict_density - target_density;
                fpreal pressure = delta * density_error;
                if(pressure < 0)
                {
                    pressure *= 0;
                    density_error *= 0;
                }
                density_error_handle.set(offset, density_error);
                fpreal last_pressure = pressure_handle.get(offset);
                fpreal new_pressure = last_pressure + pressure;
                pressure_handle.set(offset, new_pressure);
            }

            /// 5.4 accumulate pressure force
            SpikyKernel spiky(kernel_radius);
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                UT_Vector3 pressure_force = UT_Vector3(0, 0, 0);
                for (const auto &neighbor_point_id: _neighbors[i])
                {
                    long neighbor_id = neighbor_point_id.id;
                    GA_Offset neighbor_offset = gdp.pointOffset(neighbor_id);

                    UT_Vector3 neighbor_predict_position = predict_position_handle.get(neighbor_offset);
                    UT_Vector3 predict_position = predict_position_handle.get(offset);

                    fpreal neighbor_pressure = pressure_handle.get(neighbor_offset);
                    fpreal pressure = pressure_handle.get(offset);

                    UT_Vector3 dir = neighbor_predict_position - predict_position;
                    fpreal dist = (predict_position - neighbor_predict_position).length();
                    UT_Vector3 pressure_force_i = -mass_handle.get(offset) * mass_handle.get(neighbor_offset) * (pressure / (predict_position * predict_position) + neighbor_pressure / (neighbor_predict_position * neighbor_predict_position)) * spiky.gradient(dist, dir);
                    pressure_force += pressure_force_i;
                }
                pressure_force_handle.set(offset, pressure_force);
            }

            loop++;
        }

        /// correct velocity and position
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            UT_Vector3 velocity = vel_handle.get(offset);
            UT_Vector3 position = pos_handle.get(offset);

            UT_Vector3 pressure_force = pressure_force_handle.get(offset);
            UT_Vector3 external_force = external_force_handle.get(offset);

            UT_Vector3 correct_velocity = velocity + dt * (pressure_force + external_force) / mass_handle.get(offset);
            UT_Vector3 correct_position = position + dt * correct_velocity;

            /// resolve collision
            if (correct_position.x() < volumeMin.x())
            {
                correct_position.x() = volumeMin.x();
                correct_velocity.x() = 0.0;
            }
            else if (correct_position.x() > volumeMax.x())
            {
                correct_position.x() = volumeMax.x();
                correct_velocity.x() = 0.0;
            }
            if (correct_position.y() < volumeMin.y())
            {
                correct_position.y() = volumeMin.y();
                correct_velocity.y() = 0.0;
            }
            else if (correct_position.y() > volumeMax.y())
            {
                correct_position.y() = volumeMax.y();
                correct_velocity.y() = 0.0;
            }
            if (correct_position.z() < volumeMin.z())
            {
                correct_position.z() = volumeMin.z();
                correct_velocity.z() = 0.0;
            }
            else if (correct_position.z() > volumeMax.z())
            {
                correct_position.z() = volumeMax.z();
                correct_velocity.z() = 0.0;
            }

            vel_handle.set(offset, correct_velocity);
            pos_handle.set(offset, correct_position);
        }
    }
}

fpreal PCISPHSolver::_compute_delta(GU_Detail &gdp) const
{
    GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));
    GA_Offset offset = gdp.pointOffset(0); /// TODO: mass is the same for all fluid particles
    fpreal mass = mass_handle.get(offset);

    //////////////////////////////////////////////////////////////////////
    fpreal radius 			= 0.029;
    fpreal target_density 	= 1000; // water density
    fpreal target_spacing 	= radius;
    fpreal relative_kernel_radius = 1.7;
    fpreal kernel_radius 		= target_spacing * relative_kernel_radius;

    UT_Vector3 volumeMin = UT_Vector3 (-3, 0, -1);
    UT_Vector3 volumeMax = UT_Vector3 (3, 3, 1);
    fpreal h = kernel_radius;

    int iteration = 3;
    UT_Vector3 gravity(0, -9.81, 0);

    fpreal rest_density = target_density;
    fpreal dt = 0.002;
    /////////////////////////////////////////////////////////////////////


    SpikyKernel spiky(kernel_radius);
    UT_Vector3 sumGradW = UT_Vector3(0, 0, 0);
    fpreal sumGradW2 = 0;
    const fpreal supportRadius = kernel_radius;
    const fpreal diam = static_cast<fpreal>(2.0) * radius;
    const UT_Vector3 xi(0, 0, 0);
    UT_Vector3 xj = UT_Vector3(-supportRadius, -supportRadius, -supportRadius);

    while(xj.x() <= supportRadius)
    {
        while (xj.y() <= supportRadius)
        {
            while (xj.z() <= supportRadius)
            {
                fpreal dist = (xi - xj).length();
                if (dist * dist < supportRadius * supportRadius)
                {
                    UT_Vector3 dir = (xi - xj) / dist;
                    const UT_Vector3 gradW = spiky.gradient(dist, dir);
                    sumGradW += gradW;
                    sumGradW2 += gradW.dot(gradW);
                }
                xj.z() += diam;
            }
            xj.y() += diam;
            xj.z() = -supportRadius;
        }
        xj.x() += diam;
        xj.y() = -supportRadius;
        xj.z() = -supportRadius;
    }

    fpreal denom = -sumGradW.dot(sumGradW) - sumGradW2;
    fpreal beta = 2 * (mass * dt / rest_density) * (mass * dt / rest_density);
    return (std::fabs(denom) > 0) ? -1 / (beta * denom) : 0;
}


