# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

"""
    trapezoid_interp(xq, yq, xvals, yvals, zvals; compute_derivs=true) -> (z, zx, zy)

Bilinear interpolation of `z` at `(xq, yq)` inside a trapezoid defined by four
points P1=BL, P2=BR, P3=TL, P4=TR. The coordinates are provided in `xvals` and
`yvals`, and the corresponding values in `zvals` (each length 4, in the same
order). When `compute_derivs=true`, returns the partial derivatives ∂z/∂x and
∂z/∂y computed via the trapezoid parameterization; otherwise returns `NaN`s for
the derivatives.
"""
function trapezoid_interp(
        xq, yq,
        xvals, yvals, zvals;
        compute_derivs = true
    )

    # Unpack points
    x1, y1, z1 = xvals[1], yvals[1], zvals[1]
    x2, y2, z2 = xvals[2], yvals[2], zvals[2]
    x3, y3, z3 = xvals[3], yvals[3], zvals[3]
    x4, y4, z4 = xvals[4], yvals[4], zvals[4]

    # Solve for t (vertical param)
    v_left = (x3 - x1, y3 - y1)
    v_point = (xq - x1, yq - y1)
    denom = v_left[1]^2 + v_left[2]^2
    t = (v_left[1] * v_point[1] + v_left[2] * v_point[2]) / denom

    # Interpolate left and right points for this t
    xl = x1 + t * (x3 - x1)
    yl = y1 + t * (y3 - y1)
    xr = x2 + t * (x4 - x2)
    yr = y2 + t * (y4 - y2)

    # Solve for s (horizontal param)
    dx = xr - xl; dy = yr - yl
    denom2 = dx^2 + dy^2
    s = ((xq - xl) * dx + (yq - yl) * dy) / denom2

    # Interpolate z bilinearly
    z_bottom = (1 - s) * z1 + s * z2
    z_top = (1 - s) * z3 + s * z4
    z = (1 - t) * z_bottom + t * z_top

    if compute_derivs
        # Derivatives w.r.t s,t in param space
        dz_ds = (1 - t) * (z2 - z1) + t * (z4 - z3)
        dz_dt = (z_top - z_bottom)

        # Chain rule: ∂(x,y)/∂(s,t)
        dx_ds = dx; dy_ds = dy
        dx_dt = (x3 - x1) + s * ((x4 - x2) - (x3 - x1))
        dy_dt = (y3 - y1) + s * ((y4 - y2) - (y3 - y1))

        detJ = dx_ds * dy_dt - dx_dt * dy_ds
        zx = (dz_ds * dy_dt - dz_dt * dy_ds) / detJ
        zy = (-dz_ds * dx_dt + dz_dt * dx_ds) / detJ
    else
        zx, zy = NaN, NaN
    end

    return z, zx, zy
end

function find_bracketing_index(X::AbstractVector{<:Real}, xq::Real)

    n = length(X)
    # Find bracketing index
    i = searchsortedfirst(X, xq) - 1

    if i <= 0
        return 1
    elseif i >= n
        return n - 1
    end
    return i
end

function find_bracketing_indices_Ph(Pquery::Real, hquery::Real, PVector, hMatrix)

    indx_P = find_bracketing_index(PVector, Pquery)
    indx_hB = find_bracketing_index(hMatrix[indx_P, :], hquery)
    indx_hT = find_bracketing_index(hMatrix[indx_P + 1, :], hquery)

    # A tuple of 2-tuples, not a Vector of Vectors: `A[i...]` where `i` is a Vector has
    # a runtime-unknown number of indices, so it infers as `Any`. With an NTuple{2,Int}
    # the splat is statically known and the element type is recovered.
    idxs = ((indx_P, indx_hB), (indx_P, indx_hB + 1), (indx_P + 1, indx_hT), (indx_P + 1, indx_hT + 1))

    return idxs # P1=BL, P2=BR, P3=TL, P4=TR.
end

# Property tables are selected by `Val` rather than by a runtime-computed Symbol.
# `getfield(data, Symbol(string(name, :Matrix)))` is not inferable, so the old code
# returned `Tuple{Any,Any,Any}` from every property call; that `Any` propagated into
# the model equations and made them uncompilable under `juliac --trim`.
_prop_matrix(data::RefrigerantData, ::Val{:Density}) = data.DensityMatrix
_prop_matrix(data::RefrigerantData, ::Val{:Temperature}) = data.TemperatureMatrix
_prop_matrix(::RefrigerantData, ::Val) = nothing

_prop_vector(data::RefrigerantData, ::Val{:BubbleEnthalpy}) = data.BubbleEnthalpyVector
_prop_vector(data::RefrigerantData, ::Val{:DewEnthalpy}) = data.DewEnthalpyVector
_prop_vector(::RefrigerantData, ::Val) = nothing

function interp2D_property_Ph(P, h, ::Val{PROP}, data::RefrigerantData; compute_derivs = true, verbose = true) where {PROP}
    if !checkbounds_P(P, data) || !checkbounds_h(h, data) || isnan(P) || isnan(h)
        return (NaN, NaN, NaN)
    end

    propMat = _prop_matrix(data, Val(PROP))
    if propMat === nothing
        verbose && @error "Unsupported property name."
        return (NaN, NaN, NaN)
    end
    idxs = find_bracketing_indices_Ph(P, h, data.PVector, data.hMatrix)
    PVals = [data.PVector[i[1]] for i in idxs]
    hVals = [data.hMatrix[i...] for i in idxs]
    propVals = [propMat[i...] for i in idxs]

    return trapezoid_interp(P, h, PVals, hVals, propVals; compute_derivs = compute_derivs)
end

# Symbol-taking method kept for API compatibility.
interp2D_property_Ph(P, h, property_name::Symbol, data::RefrigerantData; kwargs...) =
    interp2D_property_Ph(P, h, Val(property_name), data; kwargs...)

function interp1D_property_P(P, ::Val{PROP}, data::RefrigerantData; verbose = true) where {PROP}
    if !checkbounds_P(P, data) || isnan(P)
        return (NaN, NaN)
    end

    propVec = _prop_vector(data, Val(PROP))
    if propVec === nothing
        verbose && @error "Unsupported property name."
        return (NaN, NaN)
    end

    PVector = data.PVector
    i = find_bracketing_index(PVector, P)
    t = (P - PVector[i]) / (PVector[i + 1] - PVector[i])
    return (1 - t) * propVec[i] + t * propVec[i + 1], (propVec[i + 1] - propVec[i]) / (PVector[i + 1] - PVector[i])
end

function checkbounds_P(P, data::RefrigerantData)
    if P < data.PMin || P > data.PMax
        VERBOSE && @warn("Input P=$P must be in the interval [data.PMin, data.PMax]. Returning NaNs.")
        return false
    else
        return true
    end
end

function checkbounds_h(h, data::RefrigerantData)
    if h < data.hMin || h > data.hMax
        VERBOSE && @warn("Input h=$h must be in the interval [data.hMin, data.hMax]. Returning NaNs.")
        return false
    else
        return true
    end
end

# Symbol-taking method kept for API compatibility.
interp1D_property_P(P, property_name::Symbol, data::RefrigerantData; kwargs...) =
    interp1D_property_P(P, Val(property_name), data; kwargs...)
