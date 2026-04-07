# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

using Pkg.Artifacts
using JLD2

# Custom show method for better display of RefrigerantData
function Base.show(io::IO, ::MIME"text/plain", x::RefrigerantData)
    return print(io, "$(x.name) Data")
end

function refrigerant_data_dir()
    return artifact"refrigerant-data"
end

function refrigerant_file_path(name::AbstractString)
    path = joinpath(refrigerant_data_dir(), "$(name).jld2")
    isfile(path) || error("Unknown refrigerant \"$name\". No file at $path")
    return path
end

"""
    available_refrigerants() -> Vector{String}

Return a sorted list of available refrigerant names based on the .jld2 files present in the data artifacts.
"""
function available_refrigerants()
    files = filter(f -> endswith(f, ".jld2"), readdir(refrigerant_data_dir()))
    return sort(first.(splitext.(files)))
end

"""
    has_refrigerant(name::AbstractString) -> Bool

Check if a refrigerant dataset with the given name exists in the data artifacts. Returns `true` if a corresponding .jld2 file is found, `false` otherwise.
"""
function has_refrigerant(name::AbstractString)
    return isfile(joinpath(refrigerant_data_dir(), "$(name).jld2"))
end

"""
    load_refrigerant(ref_name::AbstractString) -> RefrigerantData

Load a refrigerant dataset from the package `data/` directory using the
provided `ref_name` (expects a `((ref_name).jld2` file). The returned value is
validated with `validate_data` before being returned. Raises an error if the
file cannot be found or the stored data is invalid.
"""
function load_refrigerant(ref_name::AbstractString)
    if !has_refrigerant(ref_name)
        error("Refrigerant \"$ref_name\" not found. Available refrigerants: $(available_refrigerants())")
    end
    data_path = refrigerant_file_path(ref_name)
    @info "Loading data from $data_path"
    data = RefrigerantData(load(data_path)["data"])
    validate_data(data)
    return data
end

"""
    validate_data(data::RefrigerantData)

Validate the internal consistency of a `RefrigerantData` instance. Checks include
monotonic pressure grids, bounds consistency, matrix/vector sizes, and
dimensions matching between the property tables and the P–h grid. Raises an
`AssertionError` if any condition fails and logs a success message otherwise.
"""
function validate_data(data::RefrigerantData)
    n = length(data.PVector)
    @assert n ≥ 2 "PVector must have at least 2 entries."
    @assert size(data.hMatrix, 1) == n "hMatrix must have same number of rows as PVector."
    @assert data.PMin < data.PMax "PMin must be less than PMax."
    @assert data.hMin < data.hMax "hMin must be less than hMax."
    @assert all(data.PVector .>= data.PMin) "All PVector values must be ≥ PMin."
    @assert all(data.PVector .<= data.PMax) "All PVector values must be ≤ PMax."
    @assert all(data.hMatrix .>= data.hMin) "All hMatrix values must be ≥ hMin."
    @assert all(data.hMatrix .<= data.hMax) "All hMatrix values must be ≤ hMax."
    @assert all(diff(data.PVector) .> 0) "PVector must contain unique P values in ascending order."
    @assert all(diff(data.hMatrix, dims = 2) .> 0) "Rows of hMatrix must contain unique h values in ascending order."
    @assert size(data.DensityMatrix) == size(data.hMatrix) "DensityMatrix must have same size as hMatrix."
    @assert size(data.TemperatureMatrix) == size(data.hMatrix) "TemperatureMatrix must have same size as hMatrix."
    @assert length(data.BubbleEnthalpyVector) == n "BubbleEnthalpyVector must have same length as PVector."
    @assert length(data.DewEnthalpyVector) == n "DewEnthalpyVector must have same length as PVector."
    return @info "$(data.name) data successfully validated."
end
