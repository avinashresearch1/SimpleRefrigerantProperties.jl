# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

# Field types are deliberately left unconstrained. Pinning them to abstract types
# (AbstractFloat / AbstractVector{<:AbstractFloat} / AbstractMatrix{<:AbstractFloat})
# made every field read type-unstable, which propagated `Any` through every property
# call and made models using them uncompilable under `juliac --trim`. Naming the
# fields only lets the loader's concrete element type (Float32 for the shipped
# tables) flow through unchanged, so instances are concrete and inference succeeds.
const RefrigerantData = NamedTuple{
    (
        :name,
        :PMin,
        :PMax,
        :PCrit,
        :hMin,
        :hMax,
        :PVector,
        :hMatrix,
        :DensityMatrix,
        :TemperatureMatrix,
        :BubbleEnthalpyVector,
        :DewEnthalpyVector,
    ),
}
