using LinearAlgebra

using Base: +, -, *, /
using Plots

"""A 2-dimensional vector Cartesian coordinates."""
mutable struct Point2{T <: AbstractFloat}
	x::T
	y::T
end

# 2D vector algebra
Base.:+(p1::Point2, p2::Point2) = Point2(p1.x + p2.x, p1.y + p2.y)
Base.:-(p::Point2) = Point2(-p.x, -p.y)
Base.:-(p1::Point2, p2::Point2) = Point2(p1.x - p2.x, p1.y - p2.y)
Base.:*(p::Point2, x::T) where {T <: Real} = Point2(p.x * x, p.y * x)
Base.:/(p::Point2, x::T) where {T <: Real} = Point2(p.x / x, p.y / x)
LinearAlgebra.norm(p::Point2) = p.x^2 + p.y^2

a = Point2(2.0, 3.0)
b = Point2(1.0, 3.0)

"""A segment between two points."""
struct Segment{T <: AbstractFloat}
	start::Point2{T}
	finish::Point2{T} # end is a reserved keyword
end

s = Segment(a, b)

mutable struct Turtle{T <: AbstractFloat}
	pos::Point2{T}
	dir::Point2{T}
end

t = Turtle(Point2(0.0, 0.0), Point2(1.0, 0.0))

"""Rotate the turtle's walking direction by the given angle. Modifies in-place."""
function rotate(turtle::Turtle, angle::T) where {T <: Real}
	turtle.dir.x = turtle.dir.x * cos(angle) - turtle.dir.y * sin(angle)
	turtle.dir.y = turtle.dir.x * sin(angle) + turtle.dir.y * cos(angle)
	return nothing
end

"""Rotate a point about the origin by a given angle (in radians)."""
function rotate(p::Point2, angle::T) where {T <: Real}
	new_x = p.x * cos(angle) - p.y * sin(angle)
	new_y = p.x * sin(angle) + p.y * cos(angle)
	return Point2(new_x, new_y)
end

Plots.plot(points::Vector{<:Point2}) = plot([(p.x, p.y) for p in points])


"""
Draw the n-th iteration of a fractal by modifying each segment of the starting shape according
to `segment_function`. `segment_function` must take have two `Point2` argument and return a
`Vector{Point2}`.
"""
function draw_fractal(
	starting_shape::Vector{<:Point2},
	segment_function::Function;
	iterations::Unsigned=1,
)::Vector{Point2}
	if iterations == 0
		return starting_shape
	end

	# Each iteration requires the previous one to continue
	# points contains the Point2 coordinates of the current iteration
	points = starting_shape
	for _ in 1:iterations
		next_order = []

		# Each iteration consists of splitting each segment further
		for idx in eachindex(points)
			component = (
				idx == 1
				? segment_function(points[end], points[1])
				: segment_function(points[idx-1], points[idx])
			)
			append!(next_order, component)
		end

		points = next_order
	end

	return points
end
