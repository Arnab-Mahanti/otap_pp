#pragma once

#include "OTypes.h"
#include "Geometry.h"
#include "Ambient.h"
#include "Trajectory.h"
#include "Problem.h"
#include "MaterialLayer.h"
#include "Table.h"
#include "BC.h"

#include <memory>
#include <cassert>
#include <cmath>

// FIXME: Make OTAP::State thread safe

namespace OTAP
{
    class State
    {
    private:
        static inline State* state = nullptr;

        State() = default;
        ~State() = default;
        State(const State &) = delete;
        State &operator=(State &) = delete;

        // friend class cereal::access;

    public:
        // static inline auto materialManager = MaterialManager("I:/Arnab/Projects/MISC/OTAP++/docs/material_databse.yaml");
        std::unordered_map<std::string, std::shared_ptr<Geometry>> geometries;
        std::vector<Trajectory> trajectories;
        BCArray bcs;
        std::vector<std::shared_ptr<Material>> localMaterials;
        std::vector<std::shared_ptr<LayerStack>> layerstacks;

        std::vector<std::shared_ptr<ResponseProblem>> solutions;

        static State &GetInstance()
        {
            if (!state)
                state= new State();
            return *state;
        }

        // template <typename Ar>
        // void serialize(Ar &ar)
        // {
        //     ar(cereal::make_nvp("Geometries", geometries),
        //        cereal::make_nvp("Trajectories", trajectories),
        //        cereal::make_nvp("Boundary Conditions", bcs),
        //        cereal::make_nvp("Local Materials", localMaterials),
        //        cereal::make_nvp("Layerstacks", layerstacks),
        //        cereal::make_nvp("Solutions", solutions));
        //     // FIXME: State??
        // }

        // bool save(const std::string &filepath) const
        // {
        //     std::ofstream stream(filepath);
        //     return save(stream);
        // }

        // bool save(std::ostream &stream) const
        // {
        //     cereal::JSONOutputArchive ar(stream);
        //     ar(cereal::make_nvp("state", state));
        //     return true;
        // }

        // bool load(const std::string &filepath)
        // {
        //     std::ifstream stream(filepath);
        //     return load(stream);
        // }

        // bool load(std::istream &stream)
        // {
        //     cereal::JSONInputArchive ar(stream);
        //     ar(state);
        //     return true;
        // }
    };
    using Engine = State;
} // namespace OTAP
