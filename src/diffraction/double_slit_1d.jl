using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
y  = range(-1cm, 1cm, length=1000) # values where the irradiance is sampled on the screen
λ1 = 400nm
λ2 = 500nm # the wavelength of the wave
λ3 = 600nm
λ4 = 700nm
d  = 0.42mm # slit width
a  = 1mm # slit distance

wave_400 = Wave1D(9V / m, 0m, λ1, 0, c)
wave_500 = Wave1D(9V / m, 0m, λ2, 0, c)
wave_600 = Wave1D(9V / m, 0m, λ3, 0, c)
wave_700 = Wave1D(9V / m, 0m, λ4, 0, c)

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

I1 = irradiance_double_slit.(y, wave_400, d, a, s) .|> ustrip
I2 = irradiance_double_slit.(y, wave_500, d, a, s) .|> ustrip
I3 = irradiance_double_slit.(y, wave_600, d, a, s) .|> ustrip
I4 = irradiance_double_slit.(y, wave_700, d, a, s) .|> ustrip

plot(
	ustrip.(y),
	[I1, I2, I3, I4],
	title="Intensity profile of double slit diffraction",
	xlabel="Distance from center of slits",
	ylabel="Intensity",
	label=reduce(hcat, ["$(wl)nm" for wl in ["400", "500", "600", "700"]]),
	size=(1400, 1200),
	left_margin=(6.0, :mm),
)
