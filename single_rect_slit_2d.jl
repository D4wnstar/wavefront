using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
x  = range(-20cm, 20cm, length=1000) # values where the irradiance is sampled on the screen
y  = range(-20cm, 20cm, length=1000) # values where the irradiance is sampled on the screen
λ1 = 670nm
λ2 = 532nm # the wavelength of the wave
λ3 = 450nm
w  = 0.12mm # slit width
h  = 0.12mm # slit height

wave_green = Wave1D(9V / m, 0m, λ1, 0, c)
wave_red   = Wave1D(9V / m, 0m, λ2, 0, c)
wave_blue  = Wave1D(9V / m, 0m, λ3, 0, c)

"This assumes Fraunhofer's far field approximation where `distance` ≫ `slit_width`"
function irradiance_rectangular_slit(
	x::Quantity,
    y::Quantity,
	wave::AbstractWave,
	slit_width::Quantity,
    slit_height::Quantity,
	distance::Quantity,
)
    k = abs(wavenumber_ang(wave))
    α = k * 0.5slit_width * x / distance
	β = k * 0.5slit_height * y / distance
    I = irradiance(wave) * (sin(α) / α)^2 * (sin(β) / β)^2
	return I / exp(-ustrip(15abs(x)) - ustrip(15abs(y))) # dampens the center peak to show the tails better
end

I1 = QuantityArray(zeros(length(x), length(y)), u"kg*s^-3")
for (i, xi) in enumerate(x)
    for (j, yi) in enumerate(y)
        I1[i,j] = irradiance_rectangular_slit(xi, yi, wave_green, w, h, s)
    end
end

heatmap(
    ustrip.(x),
	ustrip.(y),
	ustrip.(I1),
	title="Intensity profile of rectangular slit diffraction",
	xlabel="x distance from center of slit",
	ylabel="y distance from center of slit",
    zlabel="Intensity",
	label="670nm",
    size=(700, 600)
)