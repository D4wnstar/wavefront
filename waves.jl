using LinearAlgebra
using DynamicQuantities.Constants: c, eps_0

# ϵ0 = 8.854e-12#F/m, permittivity of free space
# c  = 299_792_458#m/s, speed of light

abstract type AbstractWave end

mutable struct Wave <: AbstractWave
	field::Vector{Real}
	position::Vector{Real}
	wavelength::Real
	phase_shift::Real
	velocity::Real
end

mutable struct Wave1D <: AbstractWave
	field::Quantity
	position::Quantity
	wavelength::Quantity
	phase_shift::Quantity
	velocity::Quantity
end

Base.broadcastable(w::AbstractWave) = Ref(w)

wavenumber(wave::AbstractWave) = inv(wave.wavelength)
frequency(wave::AbstractWave) = wave.velocity / wave.wavelength
wavenumber_ang(wave::AbstractWave) = 2π * wavenumber(wave)
frequency_ang(wave::AbstractWave) = 2π * frequency(wave)

function electric_field(wave::AbstractWave, t::Real, r::AbstractVector{<:Real})
	wave.field *
	cos(frequency_ang(wave) * t - wavenumber_ang(wave) ⋅ (r - wave.position) + wave.phase_shift)
end

# The irradiance/intensity is measured using the (co)sine average in time
# so that <E₀>ₜ = 1/2 * E₀²
irradiance(wave::AbstractWave) = c * eps_0 * 0.5wave.field^2