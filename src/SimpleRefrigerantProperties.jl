# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

module SimpleRefrigerantProperties

global VERBOSE = true

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
