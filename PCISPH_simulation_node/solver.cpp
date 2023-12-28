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

fpreal PCISPHSolver::_halfH = 0.0;
fpreal PCISPHSolver::_kernelValueCoeff = 0.0;
fpreal PCISPHSolver::_kernelGradientCoeffA = 0.0;
fpreal PCISPHSolver::_kernelGradientCoeffB = 0.0;

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

        GA_RWAttributeRef mass_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS), 1, GA_Defaults(0));

        GA_RWAttributeRef vel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, gdp.getStdAttributeName(GEO_ATTRIBUTE_VELOCITY),
                                                      3, GA_Defaults(0));

        GA_RWAttributeRef density_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "density", 1, GA_Defaults(0));
        GA_RWAttributeRef pressure_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pressure", 1, GA_Defaults(0));

        GA_RWAttributeRef predict_position_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pP", 3, GA_Defaults(0));
        GA_RWAttributeRef predict_velocity_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pV", 3, GA_Defaults(0));
        GA_RWAttributeRef predict_density_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pDensity", 1, GA_Defaults(0));

        GA_RWAttributeRef density_error_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "Derr", 1, GA_Defaults(0));

        GA_RWAttributeRef external_accel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "eA", 3, GA_Defaults(0));
        GA_RWAttributeRef pressure_accel_ref = gdp.addFloatTuple(GA_ATTRIB_POINT, "pA", 3, GA_Defaults(0));

        GA_RWAttributeRef neighbor_num_ref = gdp.addIntTuple(GA_ATTRIB_POINT, "nNum", 1, GA_Defaults(0));

        //precomputeKernelCoefficients();

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
    }
}

/*void PCISPHSolver::precomputeKernelCoefficients() const {
    const fpreal PI = 3.14159265359;
    const fpreal h = 0.2;

    _halfH = h/2.0;	// In Monaghan2005, h=half of smoothing radius

    // Precompute value coefficient (Identical for part A and B)
    _kernelValueCoeff = 1.0 / (4.0*PI*pow(_halfH,3));

    // Precompute gradient coefficients
    _kernelGradientCoeffA = 3.0 / (4.0*PI*pow(_halfH,4));
    _kernelGradientCoeffB = -3.0 / (4.0*PI*pow(_halfH,4));
}*/

namespace
{
    // This is much faster than calling pow(val, exponent)
    inline double pow2(double val) { return val*val; }
    inline double pow3(double val) { return val*val*val; }
    inline double pow7(double val) { return val*val*val*val*val*val*val; }
}

void PCISPHSolver::solve(SIM_Object &obj, const SIM_Time &dt) const {

    //////////////////////////////////////////////////////////////////////

    UT_Vector3 volumeMin = UT_Vector3 (-0.5, 0, -0.5);
    UT_Vector3 volumeMax = UT_Vector3 (0.5, 1.5, 0.5);
    fpreal radius 			= 0.029;

    fpreal target_spacing 	= radius;
    fpreal relative_kernel_radius = 1.7;
    fpreal kernel_radius 		= target_spacing * relative_kernel_radius;
    fpreal h = kernel_radius;
    fpreal target_density 	= 1000.0; // water density

    fpreal _dt = 0.002;
    UT_Vector3 gravity(0, -9.81, 0);
    /////////////////////////////////////////////////////////////////////

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
        GA_RWHandleF density_handle = gdp.findPointAttribute("density");
        GA_RWHandleF pressure_handle = gdp.findPointAttribute("pressure");

        GA_RWHandleV3 predict_position_handle = gdp.findPointAttribute("pP");
        GA_RWHandleV3 predict_velocity_handle = gdp.findPointAttribute("pV");
        GA_RWHandleF predict_density_handle = gdp.findPointAttribute("pDensity");
        GA_RWHandleF density_error_handle = gdp.findPointAttribute("Derr");

        GA_RWHandleV3 external_accel_handle = gdp.findPointAttribute("eA");
        GA_RWHandleV3 pressure_accel_handle = gdp.findPointAttribute("pA");

        GA_RWHandleI neighbor_num_handle = gdp.findPointAttribute("nNum");

        GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));

        long long particle_size = gdp.getNumPoints();
        /// 1.update Neighbors
        std::vector<std::vector<Neighbor> >	_neighbors;
        // Reserve space and initialize neighbors' data
        _neighbors.clear();
        _neighbors.resize(particle_size);
        for (long p=0; p<particle_size; ++p)
        {
            neighbor_num_handle.set(p, 0);
        }
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
            density_handle.set(offset, 0.0);
            fpreal particle_density = mass_handle(offset) * poly6(0.0);
            for (const auto & neighbor : _neighbors[i])
            {
                fpreal dist = (pos_handle.get(offset) - pos_handle.get(gdp.pointOffset(neighbor.id))).length();
                particle_density += mass_handle(offset) * poly6(dist);
            }
            density_handle.set(offset, particle_density);
        }

        /// 3.add external force
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            UT_Vector3 gravity_accel = UT_Vector3(0.0, -9.81, 0.0);

            UT_Vector3 viscosity_accel = UT_Vector3(0.0, 0.0, 0.0);
            for(const auto & neighbor : _neighbors[i])
            {
                GA_Offset neighborOffset = gdp.pointOffset(neighbor.id);
                fpreal dist = (pos_handle.get(offset) - pos_handle.get(neighborOffset)).length();
                UT_Vector3 dir = pos_handle.get(offset) - pos_handle.get(neighborOffset);
                viscosity_accel += 0.5 * (2+2) * 0.01 * (mass_handle(offset) / density_handle.get(neighborOffset) * (vel_handle.get(offset) - vel_handle.get(neighborOffset)).dot(dir)) / (dist * dist + 0.01 * radius * radius) * poly6.cubic_kernel_derivative(dist);
            }
            UT_Vector3 external_accel = gravity_accel + viscosity_accel;
            external_accel_handle.set(offset, external_accel);
        }

        /// 4.initialize pressure and pressure force
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            pressure_handle.set(offset, 0);
            pressure_accel_handle.set(offset, UT_Vector3(0, 0, 0));
            density_error_handle.set(offset, 0);
            //predict_density_handle.set(offset, density_handle.get(offset));
        }

        /// 5.prediction and correction
        int loop = 0;
        int min_loop = 3;
        bool density_error_too_large = true;
        int max_loop = 5;
        while(((loop < min_loop)||(density_error_too_large)) && (loop < max_loop))
        {
            /// 5.1 predict velocity and position
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                UT_Vector3 velocity = vel_handle.get(offset);
                UT_Vector3 position = pos_handle.get(offset);

                UT_Vector3 predict_velocity = velocity + dt * external_accel_handle.get(offset) + dt * pressure_accel_handle.get(offset);
                UT_Vector3 predict_position = position + dt * predict_velocity;

                if(predict_position.x() < volumeMin.x())
                {
                    predict_position.x() = volumeMin.x();
                    predict_velocity.x() = 0.0;
                }
                else if(predict_position.x() > volumeMax.x())
                {
                    predict_position.x() = volumeMax.x();
                    predict_velocity.x() = 0.0;
                }
                if(predict_position.y() < volumeMin.y())
                {
                    predict_position.y() = volumeMin.y();
                    predict_velocity.y() = 0.0;
                }
                else if(predict_position.y() > volumeMax.y())
                {
                    predict_position.y() = volumeMax.y();
                    predict_velocity.y() = 0.0;
                }
                if(predict_position.z() < volumeMin.z())
                {
                    predict_position.z() = volumeMin.z();
                    predict_velocity.z() = 0.0;
                }
                else if(predict_position.z() > volumeMax.z())
                {
                    predict_position.z() = volumeMax.z();
                    predict_velocity.z() = 0.0;
                }

                predict_velocity_handle.set(offset, predict_velocity);
                predict_position_handle.set(offset, predict_position);
            }

            /// 5.2 predict density
            //StdKernel poly6(kernel_radius);
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                fpreal sum = mass_handle(offset) * poly6(0);
                for (long n = 0; n < _neighbors[i].size(); ++n)
                {
                    const Neighbor &neighbor = _neighbors[i][n];
                    GA_Offset neighborOffset = gdp.pointOffset(neighbor.id);
                    fpreal predict_dist = (predict_position_handle.get(offset) - predict_position_handle.get(neighborOffset)).length();
                    sum += mass_handle(offset) * poly6(predict_dist);
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
                fpreal new_pressure = pressure_handle.get(offset);
                new_pressure += pressure;
                pressure_handle.set(offset, new_pressure);
            }

            /// 5.4 accumulate pressure force
            SpikyKernel spiky(kernel_radius);
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                for (long n = 0; n < _neighbors[i].size(); ++n)
                {
                    const Neighbor &neighbor = _neighbors[i][n];
                    GA_Offset neighborOffset = gdp.pointOffset(neighbor.id);

                    fpreal predict_dist = (predict_position_handle.get(offset) - predict_position_handle.get(neighborOffset)).length();
                    if(predict_dist > std::numeric_limits<fpreal>::epsilon() && predict_density_handle.get(neighborOffset) > std::numeric_limits<fpreal>::epsilon())
                    {
                        UT_Vector3 predict_dir = (predict_position_handle.get(neighborOffset) - predict_position_handle.get(offset))/predict_dist;
                        // Compute contribution
                        UT_Vector3 test1 = spiky.gradient(predict_dist, predict_dir);
                        fpreal test2 = mass_handle(neighborOffset);
                        fpreal test3 = pressure_handle.get(neighborOffset);
                        fpreal test4 = predict_density_handle.get(neighborOffset);
                        fpreal test5 = pressure_handle.get(offset);
                        fpreal test6 = predict_density_handle.get(offset);
                        fpreal test7 = pressure_handle.get(offset) / pow2(predict_density_handle.get(offset));
                        fpreal test8 = pressure_handle.get(neighborOffset) / pow2(predict_density_handle.get(neighborOffset));

                        UT_Vector3 contribution = mass_handle(neighborOffset) * (pressure_handle.get(offset) / pow2(predict_density_handle.get(offset)) + pressure_handle.get(neighborOffset) / pow2(predict_density_handle.get(neighborOffset))) * spiky.gradient(predict_dist, predict_dir);
                        // Add contribution
                        UT_Vector3 pressure_accel = pressure_accel_handle.get(offset);
                        pressure_accel -= contribution;
                        pressure_accel_handle.set(offset, pressure_accel);
                    }
                }
            }

            density_error_too_large = false;
            unsigned int overCount = 0;
            for(long long i = 0; i < particle_size; i++)
            {
                GA_Offset offset = gdp.pointOffset(i);
                fpreal density_error = density_error_handle.get(offset);
                if(density_error > 0.01 * target_density)
                {
                    overCount++;
                }
            }

            if(overCount > 0)
            {
                density_error_too_large = true;
            }

            loop++;
        }

        /// correct velocity and position
        for(long long i = 0; i < particle_size; i++)
        {
            GA_Offset offset = gdp.pointOffset(i);
            UT_Vector3 velocity = vel_handle.get(offset);
            UT_Vector3 position = pos_handle.get(offset);

            UT_Vector3 pressure_accel = pressure_accel_handle.get(offset);
            UT_Vector3 external_accel = external_accel_handle.get(offset);

            velocity += dt * (pressure_accel + external_accel);
            position += dt * velocity;

            /// resolve collision
            if (position.x() < volumeMin.x())
            {
                position.x() = volumeMin.x();
                velocity.x() = 0.0;
            }
            else if (position.x() > volumeMax.x())
            {
                position.x() = volumeMax.x();
                velocity.x() = 0.0;
            }
            if (position.y() < volumeMin.y())
            {
                position.y() = volumeMin.y();
                velocity.y() = 0.0;
            }
            else if (position.y() > volumeMax.y())
            {
                position.y() = volumeMax.y();
                velocity.y() = 0.0;
            }
            if (position.z() < volumeMin.z())
            {
                position.z() = volumeMin.z();
                velocity.z() = 0.0;
            }
            else if (position.z() > volumeMax.z())
            {
                position.z() = volumeMax.z();
                velocity.z() = 0.0;
            }

            vel_handle.set(offset, velocity);
            pos_handle.set(offset, position);
        }
    }
}

fpreal PCISPHSolver::_compute_delta(GU_Detail &gdp) const
{
    fpreal mass = 1.0;

    //////////////////////////////////////////////////////////////////////
    fpreal radius 			= 0.029;
    fpreal target_density 	= 1000.0; // water density
    fpreal kernel_radius 		= 0.029 * 1.7;

    UT_Vector3 volumeMin = UT_Vector3 (-0.5, 0, -0.5);
    UT_Vector3 volumeMax = UT_Vector3 (0.5, 1.5, 0.5);
    fpreal h = kernel_radius;

    fpreal rest_density = target_density;
    fpreal dt = 0.002;

    GA_RWHandleF mass_handle = gdp.findPointAttribute(gdp.getStdAttributeName(GEO_ATTRIBUTE_MASS));
    /////////////////////////////////////////////////////////////////////

    UT_Vector3 sumGradW = UT_Vector3(0, 0, 0);
    fpreal sumGradW2 = 0.0;
    const float supportRadius = kernel_radius;
    const fpreal diam = static_cast<fpreal>(2.0) * radius;
    const UT_Vector3 xi(0, 0, 0);
    UT_Vector3 xj = {-supportRadius, -supportRadius, -supportRadius};

    SpikyKernel spiky(kernel_radius);
    while(xj[0] <= supportRadius)
    {
        while (xj[1] <= supportRadius)
        {
            while (xj[2] <= supportRadius)
            {
                fpreal dist = (xi - xj).length();
                if (dist * dist < supportRadius * supportRadius)
                {
                    UT_Vector3 dir = (xj - xi) / dist;
                    const UT_Vector3 gradW = spiky.gradient(dist, dir);
                    sumGradW += gradW;
                    sumGradW2 += gradW.dot(gradW);
                }
                xj[2] += diam;
            }
            xj[1] += diam;
            xj[2] = -supportRadius;
        }
        xj[0] += diam;
        xj[1] = -supportRadius;
        xj[2] = -supportRadius;
    }

    fpreal denom = -sumGradW.dot(sumGradW) - sumGradW2;
    fpreal beta = 2 * (mass_handle.get(0) * dt / rest_density) * (mass_handle.get(0) * dt / rest_density);
    return (std::fabs(denom) > 0) ? -1 / (beta * denom) : 0;
}

/*fpreal PCISPHSolver::getKernelValue(fpreal dist) {
    fpreal q = dist/_halfH;
    if (q<1.0)
    {
        return _kernelValueCoeff * ( pow3(2.0-q)-4*pow3(1.0-q) );
    }
    else
    {
        return _kernelValueCoeff * pow3(2.0-q);
    }
}

UT_Vector3 PCISPHSolver::getKernelGradient(fpreal dist, const UT_Vector3 &xij) {
    fpreal q = dist/_halfH;
    UT_Vector3 gradient = xij;
    if (q<= 0.0)
    {
        gradient = UT_Vector3(0,0,0);
    }
    else if (q<1.0)
    {
        gradient *= _kernelGradientCoeffA * (4.0 * pow2(1.0-q) - pow2(2.0-q)) / dist;
    }
    else
    {
        gradient *= (_kernelGradientCoeffB * pow2(2.0 - q)) / dist;
    }

    return gradient;
}*/
