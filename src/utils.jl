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

    idxs = [[indx_P, indx_hB], [indx_P, indx_hB + 1], [indx_P + 1, indx_hT], [indx_P + 1, indx_hT + 1]]

    return idxs # P1=BL, P2=BR, P3=TL, P4=TR.
end

function interp2D_property_Ph(P, h, property_name, data::RefrigerantData; compute_derivs = true, verbose = true)
    if !checkbounds_P(P, data) || !checkbounds_h(h, data) || any(isnan.([P, h]))
        return (NaN, NaN, NaN)
    end

    idxs = find_bracketing_indices_Ph(P, h, data.PVector, data.hMatrix)
    PVals = [data.PVector[i[1]] for i in idxs]
    hVals = [data.hMatrix[i...] for i in idxs]

    propKey = Symbol(string(property_name, :Matrix))
    if haskey(data, propKey)
        propVals = [getfield(data, propKey)[i...] for i in idxs]
    else
        verbose && @error "Unsupported property name."
        return (NaN, NaN, NaN)
    end
    return trapezoid_interp(P, h, PVals, hVals, propVals; compute_derivs = compute_derivs)
end

function interp1D_property_P(P, property_name, data::RefrigerantData; verbose = true)
    if !checkbounds_P(P, data) || isnan(P)
        return (NaN, NaN)
    end

    propKey = Symbol(string(property_name, :Vector))
    if haskey(data, propKey)
        propVec = getfield(data, propKey)
    else
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
        VERBOSE[] && @warn("Input P=$P must be in the interval [data.PMin, data.PMax]. Returning NaNs.")
        return false
    else
        return true
    end
end

function checkbounds_h(h, data::RefrigerantData)
    if h < data.hMin || h > data.hMax
        VERBOSE[] && @warn("Input h=$h must be in the interval [data.hMin, data.hMax]. Returning NaNs.")
        return false
    else
        return true
    end
end
