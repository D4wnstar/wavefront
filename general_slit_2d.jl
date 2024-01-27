using Plots
using PlotThemes
using DynamicQuantities
using DynamicQuantities.Units: m, cm, mm, nm, V
using DynamicQuantities.Constants: c
using Meshes
using Trapz
using ProgressMeter
import WGLMakie

theme(:dark)

include("waves.jl")

s  = 4m # distance between slit and screen
x  = range(-20cm, 20cm, length=100) # values where the irradiance is sampled on the screen
y  = range(-20cm, 20cm, length=100) # values where the irradiance is sampled on the screen
λ1 = 400nm # the wavelength of the wave
λ2 = 500nm # the wavelength of the wave
λ3 = 600nm # the wavelength of the wave
λ4 = 700nm # the wavelength of the wave

wave_400 = Wave1D(9V / m, 0m, λ1, 0, c)
wave_500 = Wave1D(9V / m, 0m, λ2, 0, c)
wave_600 = Wave1D(9V / m, 0m, λ3, 0, c)
wave_700 = Wave1D(9V / m, 0m, λ4, 0, c)

waves = [wave_400, wave_500, wave_600, wave_700]
# waves = [wave_400]

# verts = [
# 	(0.0, 0.0),
# 	(0.0, 1.0e-4),
# 	(1.0e-4, 1.0e-4),
# 	(2.0e-4, 2.0e-4),
# 	(1.0e-4, 3.0e-4),
# 	(0.0, 2.0e-4),
# 	(-1.0e-4, 1.0e-4),
# ]
verts = [
	(-6e-4, -6e-4),
	(-6e-4, 6e-4),
	(6e-4, -6e-4),
	(6e-4, 6e-4),
]
slit = Ngon(verts...)

function make_grid_in_geometry(f::Function, geometry::Ngon)
	box = boundingbox(geometry)
	minx, miny = coordinates(box.min)
	maxx, maxy = coordinates(box.max)
	x = range(minx, maxx, length=100)
	y = range(miny, maxy, length=100)
	grid = Iterators.product(x, y)
	vals = map(p -> Point(ustrip.(p)) ∈ geometry ? f(p) : 0.0, grid)
	return x, y, collect(vals') # The resulting mapped array is transposed compared to the original
end

# The electric field E of a wave going through a generic slit of area A is
# E = e^iks / A * double integral of incident field * e^-i(k_x*x'+k_y*y') dx'dy'
# with k_x = ksin(θ)cos(ϕ) and k_y = ksin(θ)cos(ϕ)
# The integral is over the area of the slit
# With a sufficiently far away screen (at distance s), sin(θ) ≈ tan(θ) ≈ y/s
# and cos(θ) ≈ 1 + tan(θ) ≈ 1 + y/s. Same goes for ϕ, except it becomes x/s
# Thus we get k_x = ky/s + kyx/s^2 and k_y = kyx/s^2
# Calling yx/s^2 = β and yx/s^2 + y/s = α and incident field E0 we get the integral
# ∬E0(x',y')e^[-ik(αx'+βy')] dx'dy'
# which is the (2D) Fourier transform of E0, so Ê0(αk,βk)
# Thus the electric field for over point (x,y) of the screen is
# E(x,y) = e^iks / A * Ê0(αk, βk)
# Partly from Wikipedia (https://en.wikipedia.org/wiki/Diffraction#General_aperture)
# partly derived myself

"This assumes Fraunhofer's far field approximation"
function irradiance_general_slit(
	x::Quantity,
	y::Quantity,
	wave::AbstractWave,
	slit_ngon::Ngon,
	s::Quantity,
)
	k = ustrip(wavenumber_ang(wave))
	A = area(slit_ngon)
	E_0 = ustrip(wave.field)
	I_0 = ustrip(irradiance(wave))
	box = boundingbox(slit)
	minx, miny = coordinates(box.min)
	maxx, maxy = coordinates(box.max)
	wx = abs(maxx - minx)
	wy = abs(maxy - miny)
	xc = 0.5(maxx + minx)
	yc = 0.5(maxy + miny)
	β = x / s
	α = y / s
	a = wx*α*k
	b = wy*β*k
	# I = 4E_0 / (A * k^2) * sin(wx * α * k) * sin(wy * β * k) * exp(-im * k * (xc + yc))
	I = I_0 * (sin(a) / a)^2 * (sin(b) / b)^2 * exp(-im * k * (xc + yc))
	# Ex, Ey, E::Matrix{Number} = make_grid_in_geometry(
	# 	((x, y),) -> real(E_0 * exp(-im * k * (α * x + β * y))),
	# 	slit_ngon,
	# )
	# I = E_0 / A * trapz((Ex, Ey), E)
	return I# / exp(ustrip(-15abs(x) - 15abs(y))) # dampens the center peak to show the tails better
end

I = [zeros(length(x), length(y)) for _ in 1:4]
p = Progress(length(x) * length(y), desc="Computing irradiance...")
Threads.@threads for widx in eachindex(waves)
	for i in eachindex(x)
		for j in eachindex(y)
			I[widx][i, j] = irradiance_general_slit(x[i], y[i], waves[widx], slit, s) |> real
			next!(p)
		end
	end
end
finish!(p)


heatmap(
	ustrip.(x),
	ustrip.(y),
	I,
	title=reshape(
		[
			"Intensity profile of rectangular slit diffraction [$(wl)nm]" for
			wl in ["400", "500", "600", "700"]
		], (1, 4)),
	xlabel="x distance from center of slit",
	ylabel="y distance from center of slit",
	colorbar_title="Intensity",
	left_margin=(6.0, :mm),
	size=(1400, 1200),
	layout=4,
	aspectratio=1,
)
