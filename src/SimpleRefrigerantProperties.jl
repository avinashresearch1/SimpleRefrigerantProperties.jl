# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

module SimpleRefrigerantProperties

# Toggle for emitting @warn messages on out-of-bounds inputs. Stored in a
# `const Ref` so it can be mutated (e.g. silenced during tests) without
# triggering Julia's non-const global / cross-module binding warnings.
const VERBOSE = Ref(true)

include("types.jl")
include("data_utils.jl")
include("utils.jl")
include("properties.jl")

export RefrigerantData

export available_refrigerants,
    has_refrigerant,
    load_refrigerant,
    validate_data

export DerDensity_Ph,
    DerTemperature_Ph,
    BubbleEnthalpy_P,
    DewEnthalpy_P
end
