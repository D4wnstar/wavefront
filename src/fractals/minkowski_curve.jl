using Meshes
using Rotations
using LinearAlgebra
using Plots
using PlotThemes

theme(:dracula)

include("common.jl")

function segment_minkowski_component(A::Point2, Z::Point2)::Vector{Point2}
	d = norm(Z - A)
	dir = (Z - A) / d
	step = d / 4

	B = A + dir * step
	C = B + rotate(dir, π / 2) * step
	D = C + dir * step
	E = D + rotate(dir, -π / 2) * step
	F = E + rotate(dir, -π / 2) * step
	G = F + dir * step
	H = G + rotate(dir, π / 2) * step

	return [A, B, C, D, E, F, G, H, Z]
end

A = Point2(0.0, 0.0)
B = Point2(4.0, 0.0)

verts = [A, B]
graphs = [plot([A, B])]
titles = ["Starting shape"]
for i in 1:5
	points = draw_fractal(verts, segment_minkowski_component, iterations=unsigned(i))
	push!(graphs, plot(points))
	titles = [titles... "Iteration $i"]
end

l = @layout [a b; c d; e f]
plot(
	graphs...,
	title=titles,
	layout=l,
	size=(1200, 1500),
	legend=false,
	plot_title="Minkowski Curve",
)
