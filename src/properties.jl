# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

"""
    DerDensity_Ph(P, h, data::RefrigerantData) -> (ρ, ∂ρ∂P, ∂ρ∂h)

Interpolate density at pressure `P` and specific enthalpy `h` from the refrigerant
property tables in `data`. Returns the density along with partial derivatives
∂ρ/∂P and ∂ρ/∂h computed from the trapezoidal interpolation surface. If inputs
are out of bounds, the returned values are `NaN`.
"""
function DerDensity_Ph(P, h, data::RefrigerantData)
    return interp2D_property_Ph(P, h, Val(:Density), data; compute_derivs = true)
end

"""
    DerTemperature_Ph(P, h, data::RefrigerantData) -> (T, ∂T∂P, ∂T∂h)

Interpolate temperature at pressure `P` and specific enthalpy `h` from the
refrigerant property tables in `data`. Returns temperature along with partial
derivatives ∂T/∂P and ∂T/∂h from the trapezoidal interpolation surface. If inputs
are out of bounds, the returned values are `NaN`.
"""
function DerTemperature_Ph(P, h, data::RefrigerantData)
    return interp2D_property_Ph(P, h, Val(:Temperature), data; compute_derivs = true)
end

"""
    BubbleEnthalpy_P(P, data::RefrigerantData) -> hB

Interpolate bubble-point enthalpy at pressure `P` using the 1D property vector
from `data`. Returns `NaN` if `P` is out of bounds.
"""
function BubbleEnthalpy_P(P, data::RefrigerantData)
    return interp1D_property_P(P, Val(:BubbleEnthalpy), data)[1]
end

"""
    DewEnthalpy_P(P, data::RefrigerantData) -> hD

Interpolate dew-point enthalpy at pressure `P` using the 1D property vector
from `data`. Returns `NaN` if `P` is out of bounds.
"""
function DewEnthalpy_P(P, data::RefrigerantData)
    return interp1D_property_P(P, Val(:DewEnthalpy), data)[1]
end
