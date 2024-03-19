using Meshes
using Rotations
using LinearAlgebra
using Plots
using PlotThemes
import WGLMakie

theme(:dark)

function segment_minkowski_component(A::Point2, Z::Point2)::Vector{Point2}
    d = norm(Z - A)
    dir = (Z - A) / d
    step = d / 4
    
    B = A + dir * step
    C = B + Rotate(Angle2d(π / 2))(dir) * step
    D = C + dir * step
    E = D + Rotate(Angle2d(-π / 2))(dir) * step
    F = E + Rotate(Angle2d(-π / 2))(dir) * step
    G = F + dir * step
    H = G + Rotate(Angle2d(π / 2))(dir) * step
    
    return [A, B, C, D, E, F, G, H, Z]
end

function create_minkowski_curve(A::Point2, Z::Point2; iterations=1)::Vector{Point2}
    if iterations == 0
        return [A, Z]
    end
    
    points = []
    for order in 1:iterations
        next_order = []
		iter = order == 1 ? [A, Z] : points

        for (idx, _) in enumerate(iter)
			component = (idx == 1
						 ? segment_minkowski_component(iter[end], iter[1])
						 : segment_minkowski_component(iter[idx-1], iter[idx]))
			append!(next_order, component)
		end

		points = next_order
    end

    return points
end

A = Point2(0., 0.)
Z = Point2(4., 0.)

l = @layout [a b; c d; e f]
graphs = []
titles = []
for i in 1:6
    points = create_minkowski_curve(A, Z, iterations = i)
    tuples = [(coordinates(p)[1], coordinates(p)[2]) for p in points]
    push!(graphs, plot(tuples))
    titles = [titles... "Iteration $i"]
end
plot(
    graphs...,
    title=titles,
    layout = l,
    size=(1200, 1500)
)

# viz(Rope(create_minkowski_curve(A, Z, iterations=2)...))