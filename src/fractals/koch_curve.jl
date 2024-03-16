using Meshes
using Rotations: Angle2d
using LinearAlgebra
using WGLMakie: WGLMakie
using Plots
using PlotThemes

theme(:dracula)

# d = 1
# direction = Vec(d, 0)
# start = Point(0, 0)

# forward(point::Point2, direction::Vec2) = point + direction
# rotate(direction::Vec2, angle::Angle2d) = Rotate(angle)(direction)

# function draw_koch_component(A::Point2, direction::Vec2)
# 	B = forward(A, direction)
# 	direction = rotate(direction, Angle2d(π / 3))
# 	C = forward(B, direction)
# 	direction = rotate(direction, Angle2d(-2π / 3))
# 	D = forward(C, direction)
# 	direction = rotate(direction, Angle2d(π / 3))
# 	E = forward(D, direction)

# 	return [A, B, C, D, E], direction
# end

function segment_koch_component(A::Point2, Z::Point2)::Vector{Point2}
	unnormed = (Z - A)
	dist = round(norm(unnormed); digits = 6)
	direction = unnormed / dist
	B = A + direction * dist / 3
	D = A + direction * 2dist / 3
	E = Z
	direction = Rotate(Angle2d(-π / 3))(direction)
	C = B + direction * dist / 3

	return [A, B, C, D, E]
end

function create_koch_snowflake(triangle::Vector{Point2}; iterations = 1)::Vector{Point2}
    if length(triangle) != 3
        error("The input vector must contain exactly 3 points (a triangle)")
    end
    if iterations == 0
        return triangle
    end
    
	points = []
	for order in 1:iterations
		next_order = []
		iter = order == 1 ? triangle : points

		for (idx, _) in enumerate(iter)
			component = (idx == 1
						 ? segment_koch_component(iter[end], iter[1])
						 : segment_koch_component(iter[idx-1], iter[idx]))
			append!(next_order, component)
		end

		points = next_order
	end

	return points
end

l = @layout [a b; c d; e f]
verts = [Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(0.5, sqrt(3) / 2)]
graphs = []
titles = []
for i in 1:6
    points = create_koch_snowflake(verts, iterations = i)
    tuples = [(coordinates(p)[1], coordinates(p)[2]) for p in points]
    push!(graphs, plot(tuples))
    titles = [titles... "Iteration $i"]
end
plot(
    graphs...,
    title=titles,
    layout = l,
    size=(1000, 1500)
)