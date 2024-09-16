#pragma once
#include "OTypes.h"
#include "Table.h"

namespace OTAP
{
    // TODO: Implement axial location wise variation of BCs

    enum class BCLocation
    {
        none,
        front,
        back,
        both
    };

    struct BCBase
    {
        const BCType type;
        BCLocation location = BCLocation::none;
        BCBase(BCType ptype) : type(ptype) {}
        virtual ~BCBase() = default;
    };

    struct RadiationBC : public BCBase
    {
        Table<double> Tambient;
        RadiationBC() : BCBase(BCType::Radiation) {}
    };

    struct ConvectionBC : public BCBase
    {
        Table<double> h;
        Table<double> tg;
        ConvectionBC() : BCBase(BCType::Convection) {}
    };

    struct NaturalConvectionBC : public BCBase
    {
        double l;
        Table<double> tg;
        NaturalConvectionBC() : BCBase(BCType::NaturalConvection) { location = BCLocation::back; }
    };

    struct FluxBC : public BCBase
    {
        Table<double> flux;
        FluxBC() : BCBase(BCType::Flux) {}
    };

    // TODO: Implement properly with spatial variation and update with moving grid
    struct HeatGenBC : public BCBase
    {
        Table<double> qdot;
        HeatGenBC() : BCBase(BCType::HeatGeneration) {}
    };

    struct PropellantMassBC : public BCBase
    {
        Table<double> propmass;
        PropellantMassBC() : BCBase(BCType::PropellantMass) { location = BCLocation::back; }
    };

    static inline std::shared_ptr<BCBase> make_BC(BCType T)
    {
        switch (T)
        {
        case BCType::Convection:
            return safe_make_shared<ConvectionBC>();
        case BCType::Radiation:
            return safe_make_shared<RadiationBC>();
        case BCType::Flux:
            return safe_make_shared<FluxBC>();
        case BCType::NaturalConvection:
            return safe_make_shared<NaturalConvectionBC>();
        case BCType::HeatGeneration:
            return safe_make_shared<HeatGenBC>();
        case BCType::PropellantMass:
            return safe_make_shared<PropellantMassBC>();
        default:
            assert(false);
            return nullptr;
        }
    }
    using BC = std::shared_ptr<BCBase>;
    using BCArray = std::vector<BC>;

} // namespace OTAP
