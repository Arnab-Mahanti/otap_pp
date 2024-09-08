#pragma once

#include "OTypes.h"
#include "Geometry.h"
#include "Ambient.h"
#include "Trajectory.h"
#include "Problem.h"
#include "MaterialLayer.h"
#include "Table.h"

#include <memory>
#include <cassert>
#include <cmath>

// FIXME: Make OTAP::State thread safe

namespace OTAP
{
    class State
    {
    private:
        static inline State *state = nullptr;

        State() = default;
        ~State() = default;
        State(const State &) = delete;
        State &operator=(State &) = delete;

    public:
        Geometry geometry;
        Trajectory trajectory;

        static State &GetInstance()
        {
            if (!state)
                state = new State();
            return *state;
        }
    };
    using Engine = State;
} // namespace OTAP
