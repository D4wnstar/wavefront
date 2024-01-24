using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
y  = range(-1cm, 1cm, length=1000) # values where the irradiance is sampled on the screen
λ1 = 670nm
λ2 = 532nm # the wavelength of the wave
λ3 = 450nm
d  = 0.42mm # slit width
a  = 1mm # slit distance

wave_green = Wave1D(9V / m, 0m, λ1, 0, c)
wave_red   = Wave1D(9V / m, 0m, λ2, 0, c)
wave_blue  = Wave1D(9V / m, 0m, λ3, 0, c)

"""
This assumes Fraunhofer's far field approximation where `distance` ≫ `slit_width`.
Also, the slits are both the same width. `slit_distance` refers to the length of
space in between the two slits.
"""
function irradiance_double_slit(
	x::Quantity,
	wave::AbstractWave,
	slit_width::Quantity,
	slit_distance::Quantity,
	distance::Quantity,
)
	k = abs(wavenumber_ang(wave))
	β = k * 0.5slit_width * x / distance
	d = (2slit_width + slit_distance) / 2
	return irradiance(wave) * (sin(β) / β)^2 * cos(k * slit_width * x / distance)^2
end

I1 = irradiance_double_slit.(y, wave_green, d, a, s)
I2 = irradiance_double_slit.(y, wave_red, d, a, s)
I3 = irradiance_double_slit.(y, wave_blue, d, a, s)
plot(
	ustrip.(y),
	ustrip.(I1),
	title="Intensity profile of double slit diffraction",
	xlabel="Distance from center of slits",
	ylabel="Intensity",
	label="670nm",
)
plot!(ustrip.(y), ustrip.(I2), label="532nm")
plot!(ustrip.(y), ustrip.(I3), label="450nm")
