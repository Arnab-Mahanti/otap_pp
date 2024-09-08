#ifndef _TRAJECTORY_H
#define _TRAJECTORY_H

// #include "otap.h"
#include "OTypes.h"
#include <vector>
#include <variant>
#include <optional>
#include <memory>
#include <cassert>

#include "Ambient.h"
#include "Fluid.h"

namespace OTAP
{
    class TrajectoryBase
    {
    public:
        virtual ~TrajectoryBase() = default;
        virtual double GetFreeStream() const { assert(false); };
        virtual TimeSeries GetPinf(TimePoints) const { assert(false); };
        virtual TimeSeries GetTinf(TimePoints) const { assert(false); };
        virtual TimeSeries GetRhoinf(TimePoints) const { assert(false); };
        virtual TimeSeries GetMinf(TimePoints) const { assert(false); };
        virtual TimeSeries GetVinf(TimePoints) const { assert(false); };
        virtual TimePoints GetTimePoints() const { assert(false); };

        // virtual double GetFreeStream() const { assert(false); };
        virtual double GetPinf(double) const { assert(false); };
        virtual double GetTinf(double) const { assert(false); };
        virtual double GetRhoinf(double) const { assert(false); };
        virtual double GetMinf(double) const { assert(false); };
        virtual double GetVinf(double) const { assert(false); };

        // virtual void helpGui(){assert(false);};
    };

    class VelocityTrajectory : public TrajectoryBase
    {
        TimePoints m_timepoints;
        TimeSeries m_altitude;
        TimeSeries m_velocity;
        TimeSeries m_alpha;

        std::shared_ptr<AmbientBase> m_ambient;
        std::shared_ptr<FluidBase> m_fluid;
        void ReadTrajectory(const std::string &filename, const bool delim_whatispace = true);

    public:
        VelocityTrajectory(const AmbientType &ambientType, const FluidType &fluidType, TimePoints &time, const TimeSeries &altitude, const TimeSeries &velocity, const TimeSeries &alpha)
            : m_timepoints(time), m_altitude(altitude), m_velocity(velocity), m_alpha(alpha), m_ambient(make_ambient(ambientType)), m_fluid(make_fluid(fluidType))
        {
        }

        VelocityTrajectory(const AmbientType &ambientType, const FluidType &fluidType, const std::string &filename);

        virtual double GetFreeStream() const override;
        virtual TimeSeries GetPinf(TimePoints) const override;
        virtual TimeSeries GetTinf(TimePoints) const override;
        virtual TimeSeries GetRhoinf(TimePoints) const override;
        virtual TimeSeries GetMinf(TimePoints) const override;
        virtual TimeSeries GetVinf(TimePoints) const override;
        virtual TimePoints GetTimePoints() const override;

        virtual double GetPinf(double) const override;
        virtual double GetTinf(double) const override;
        virtual double GetRhoinf(double) const override;
        virtual double GetMinf(double) const override;
        virtual double GetVinf(double) const override;
    };

    // NOTE: Disabled template system

    // template <TrajectoryType T, typename... Args>
    // std::unique_ptr<TrajectoryBase> make_trajectory(Args &&...args)
    // {
    //     if constexpr (T == TrajectoryType::Velocity)
    //         return std::make_unique<VelocityTrajectory>(std::forward<Args>(args)...);
    //     else if constexpr (T == TrajectoryType::Mach)
    //         return std::make_unique<MachTrajectory>(std::forward<Args>(args)...);
    //     else if constexpr (T == TrajectoryType::Flux)
    //         return std::make_unique<FluxTrajectory>(std::forward<Args>(args)...);
    //     else if constexpr (T == TrajectoryType::Wind_Tunnel)
    //         return std::make_unique<TunnelTrajectory>(std::forward<Args>(args)...);
    //     else if constexpr (T == TrajectoryType::CFD_Data)
    //         return std::make_unique<CFDTrajectory>(std::forward<Args>(args)...);
    //     else
    //         static_assert(false);
    // };

    template <typename... Args>
    std::shared_ptr<TrajectoryBase> make_trajectory(TrajectoryType T, Args &&...args)
    {
        switch (T)
        {
        case TrajectoryType::Velocity:
            return std::make_shared<VelocityTrajectory>(std::forward<Args>(args)...);
        default:
            assert(false);
            return nullptr;
        }
    }
    using Trajectory = std::shared_ptr<TrajectoryBase>;

} // namespace OTAP

#endif // _TRAJECTORY_H