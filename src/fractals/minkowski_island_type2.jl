# An "island" version of the more type 2 Minkowski curve

using LinearAlgebra
using Plots
using PlotThemes

theme(:dracula)

include("common.jl")

function minkowski_segmenter(seg::Segment)::Chain
	seg_length = norm(seg)
	direction = (seg.finish - seg.start) / seg_length
	step = seg_length / 4

	turtle = Turtle(seg.start, direction)

	# Production rule is → becomes →↺→↻→↻→→↺→↺→↻→
	A = seg.start
	B = forward(turtle, step)
	rotate(turtle, π / 2)
	C = forward(turtle, step)
	rotate(turtle, -π / 2)
	D = forward(turtle, step)
	rotate(turtle, -π / 2)
	E = forward(turtle, step)
	F = forward(turtle, step)
	rotate(turtle, π / 2)
	G = forward(turtle, step)
	rotate(turtle, π / 2)
	H = forward(turtle, step)
	rotate(turtle, -π / 2)
	I = forward(turtle, step)

	return Chain([A, B, C, D, E, F, G, H, I])
end

# Box start
A = Point2(0.0, 0.0)
B = Point2(4.0, 0.0)
C = Point2(4.0, 4.0)
D = Point2(0.0, 4.0)

# Uncomment for a diamond start
# A = Point2(0.0, 0.0)
# B = Point2(2.0, 2.0)
# C = Point2(4.0, 0.0)
# D = Point2(2.0, -2.0)


box = Polygon([A, B, C, D])
graphs = [plot(box, aspectratio=1)]
titles = ["Starting shape"]
for i in 1:5
	chain = draw_fractal(box, minkowski_segmenter, i)
	push!(graphs, plot(chain, aspectratio=1))
	titles = [titles... "Iteration $i"]
end

l = @layout [a b; c d; e f]
plot(
	graphs...,
	title=titles,
	layout=l,
	size=(1200, 1500),
	legend=false,
	plot_title="Minkowski Island (Type 2)",
)
