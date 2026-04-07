# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

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
    Tuple{
        String,
        AbstractFloat,
        AbstractFloat,
        AbstractFloat,
        AbstractFloat,
        AbstractFloat,
        AbstractVector{<:AbstractFloat},
        AbstractMatrix{<:AbstractFloat},
        AbstractMatrix{<:AbstractFloat},
        AbstractMatrix{<:AbstractFloat},
        AbstractVector{<:AbstractFloat},
        AbstractVector{<:AbstractFloat},
    },
}
