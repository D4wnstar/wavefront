using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c
using Meshes
import WGLMakie

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
x  = range(-20cm, 20cm, length=1000) # values where the irradiance is sampled on the screen
y  = range(-20cm, 20cm, length=1000) # values where the irradiance is sampled on the screen
λ1 = 400nm # the wavelength of the wave
λ2 = 500nm # the wavelength of the wave
λ3 = 600nm # the wavelength of the wave
λ4 = 700nm # the wavelength of the wave
w  = 0.12mm # slit width
h  = 0.12mm # slit height

wave_400 = Wave1D(9V / m, 0m, λ1, 0, c)
wave_500 = Wave1D(9V / m, 0m, λ2, 0, c)
wave_600 = Wave1D(9V / m, 0m, λ3, 0, c)
wave_700 = Wave1D(9V / m, 0m, λ4, 0, c)

waves = [wave_400, wave_500, wave_600, wave_700]

vertices = [
    (0.0, 0.0),
    (0.0, 1.0),
    (1.0, 1.0),
    (2.0, 2.0),
    (1.0, 3.0),
    (0.0, 2.0),
    (-1.0, 1.0),
]

slit = Ngon(vertices...)
viz(slit)

area(slit)

# Electric field of a wave going through a slit is
# Incident field / Area of slit * double integral of e^-i(k_x*x'+k_y*y') dx'dy'
# with k_x = ksin(θ)cos(ϕ) and k_y = ksin(θ)cos(ϕ)
# The incident field is assumed constant across the surface of the slit
# Partly from Wikipedia (https://en.wikipedia.org/wiki/Diffraction#General_aperture)

# "This assumes Fraunhofer's far field approximation where `distance` ≫ `slit_width` and `slit_height`"
# function irradiance_general_slit(
# 	x::Quantity,
# 	y::Quantity,
# 	wave::AbstractWave,
# 	slit_width::Quantity,
# 	slit_height::Quantity,
# 	distance::Quantity,
# )
# 	k = abs(wavenumber_ang(wave))
# 	α = k * 0.5slit_width * x / distance
# 	β = k * 0.5slit_height * y / distance
# 	I = irradiance(wave) * (sin(α) / α)^2 * (sin(β) / β)^2
# 	return I / exp(ustrip(-15abs(x) - 15abs(y))) # dampens the center peak to show the tails better
# end

# I = [zeros(length(x), length(y)) for _ in 1:4]
# for (matrix, wave) in zip(I, waves)
# 	for (i, xi) in enumerate(x)
# 		for (j, yi) in enumerate(y)
# 			matrix[i, j] = irradiance_general_slit(xi, yi, wave, w, h, s) |> ustrip
# 		end
# 	end
# end

# heatmap(
# 	ustrip.(x),
# 	ustrip.(y),
# 	I,
# 	title=reshape([
# 		"Intensity profile of rectangular slit diffraction [$(wl)nm]" for
# 		wl in ["400", "500", "600", "700"]
# 	], (1, 4)),
# 	xlabel="x distance from center of slit",
# 	ylabel="y distance from center of slit",
# 	colorbar_title="Intensity",
# 	left_margin=(6.0, :mm),
# 	size=(1400, 1200),
# 	layout=4,
# 	aspectratio=1,
# )
