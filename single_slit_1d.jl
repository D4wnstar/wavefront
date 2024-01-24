using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
y  = range(-2cm, 2cm, length=1000) # values where the irradiance is sampled on the screen
λ1 = 670nm
λ2 = 532nm # the wavelength of the wave
λ3 = 450nm
d  = 0.42mm # slit width

wave_green = Wave1D(9V / m, 0m, λ1, 0, c)
wave_red   = Wave1D(9V / m, 0m, λ2, 0, c)
wave_blue  = Wave1D(9V / m, 0m, λ3, 0, c)

"This assumes Fraunhofer's far field approximation where `distance` ≫ `slit_width`"
function irradiance_single_slit(
	x::Quantity,
	wave::AbstractWave,
	slit_width::Quantity,
	distance::Quantity,
)
	β = abs(wavenumber_ang(wave)) * 0.5slit_width * x / distance
	return irradiance(wave) * (sin(β) / β)^2
end

I1 = irradiance_single_slit.(y, wave_green, d, s)
I2 = irradiance_single_slit.(y, wave_red, d, s)
I3 = irradiance_single_slit.(y, wave_blue, d, s)
plot(
	ustrip.(y),
	ustrip.(I1),
	title="Intensity profile of single slit diffraction",
	xlabel="Distance from center of slit",
	ylabel="Intensity",
	label="670nm",
)
plot!(ustrip.(y), ustrip.(I2), label="532nm")
plot!(ustrip.(y), ustrip.(I3), label="450nm")
