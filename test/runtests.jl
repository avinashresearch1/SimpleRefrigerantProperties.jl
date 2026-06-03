# Copyright (C) 2026 Mitsubishi Electric Research Laboratories (MERL)
#
# SPDX-License-Identifier: MIT

using SimpleRefrigerantProperties
using Test

SimpleRefrigerantProperties.VERBOSE = false # suppress warnings during testing

# Define query points for testing.
Pquery = 3.0e6
hquery = 400.0e3

@testset begin
    refrigerants = available_refrigerants()
    @test !isempty(refrigerants)
    for name in refrigerants # list of supported refrigerants to test

        @info "Running tests for refrigerant: $name"
        data = load_refrigerant(name)

        # Tests for user-facing functions. These tests check for NaNs, correct values at boundaries, and consistency with the input data.
        @test !any(isnan.(DerDensity_Ph(Pquery, hquery, data)))
        @test !any(isnan.(DerTemperature_Ph(Pquery, hquery, data)))
        @test !isnan(BubbleEnthalpy_P(Pquery, data))
        @test !isnan(DewEnthalpy_P(Pquery, data))

        @test all(isnan.(DerDensity_Ph(1.01 * data.PMax, hquery, data)))
        @test all(isnan.(DerTemperature_Ph(0.99 * data.PMin, hquery, data)))
        @test all(isnan.(DerDensity_Ph(Pquery, 0.99 * data.hMin, data)))
        @test all(isnan.(DerTemperature_Ph(Pquery, 1.01 * data.hMax, data)))
        @test isnan(BubbleEnthalpy_P(1.01 * data.PMax, data))
        @test isnan(DewEnthalpy_P(0.99 * data.PMin, data))

        @test DerDensity_Ph(data.PMin, data.hMin, data)[1] ≈ data.DensityMatrix[1, 1]
        @test DerDensity_Ph(data.PMin, data.hMax, data)[1] ≈ data.DensityMatrix[1, end]
        @test DerDensity_Ph(data.PMax, data.hMax, data)[1] ≈ data.DensityMatrix[end, end]
        @test DerDensity_Ph(data.PMax, data.hMin, data)[1] ≈ data.DensityMatrix[end, 1]

        @test DerTemperature_Ph(data.PMin, data.hMin, data)[1] ≈ data.TemperatureMatrix[1, 1]
        @test DerTemperature_Ph(data.PMin, data.hMax, data)[1] ≈ data.TemperatureMatrix[1, end]
        @test DerTemperature_Ph(data.PMax, data.hMax, data)[1] ≈ data.TemperatureMatrix[end, end]
        @test DerTemperature_Ph(data.PMax, data.hMin, data)[1] ≈ data.TemperatureMatrix[end, 1]

        @test BubbleEnthalpy_P(data.PMin, data) ≈ data.BubbleEnthalpyVector[1]
        @test BubbleEnthalpy_P(data.PMax, data) ≈ data.BubbleEnthalpyVector[end]
        @test DewEnthalpy_P(data.PMin, data) ≈ data.DewEnthalpyVector[1]
        @test DewEnthalpy_P(data.PMax, data) ≈ data.DewEnthalpyVector[end]

        # Tests for internal functions and code coverage only. Functions in following tests are not intended to be user-facing.
        @test !isnan(SimpleRefrigerantProperties.interp2D_property_Ph(Pquery, hquery, :Density, data; compute_derivs = false)[1])
        @test all(isnan.(SimpleRefrigerantProperties.interp2D_property_Ph(Pquery, hquery, :Density, data; compute_derivs = false)[2:end]))
        @test SimpleRefrigerantProperties.find_bracketing_index(data.PVector, 1.01 * data.PMax) == length(data.PVector) - 1
        @test all(isnan.(SimpleRefrigerantProperties.interp2D_property_Ph(Pquery, hquery, :InvalidPropertyName, data; verbose = false)))
        @test all(isnan.(SimpleRefrigerantProperties.interp1D_property_P(Pquery, :InvalidPropertyName, data; verbose = false)))

        println(" ")
    end
end
