using LinearAlgebra
using Plots
using PlotThemes

using Base: +, -, *, /

theme(:dracula)

include("common.jl")

function segment_koch_component(A::Point2, Z::Point2; antiflake=false)::Vector{Point2}
	spin = antiflake ? -1 : 1
	delta = Z - A
	dist = norm(delta)
	direction = delta / dist
	B = A + direction * dist / 3
	D = A + direction * 2dist / 3
	E = Z
	direction = rotate(direction, spin * π / 3)
	C = B + direction * dist / 3

	return [A, B, C, D, E]
end

A = Point2(0.0, 0.0)
B = Point2(0.5, sqrt(3) / 2)
C = Point2(1.0, 0.0)

verts = [A, B, C]
graphs = [plot([A, B, C, A])]
titles = ["Starting shape"]
for i in 1:5
	points = draw_fractal(verts, segment_koch_component, iterations=unsigned(i))
	push!(graphs, plot(points))
	titles = [titles... "Iteration $i"]
end

l = @layout [a b; c d; e f]
plot(
	graphs...,
	title=titles,
	layout=l,
	size=(1000, 1500),
	legend=false,
	plot_title="Koch Snowflake",
)
