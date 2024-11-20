using LinearAlgebra, CSV, DataFrames, Plots, MLDatasets

function find_nearest(x::Vector, W::Array)
	i_star = 1
	j_star = 1
	shortest_dist = Inf
	for i in axes(W, 1)
		for j in axes(W, 2)
			curr_dist = norm(x - W[i, j])
			if curr_dist < shortest_dist
				shortest_dist = curr_dist
				i_star = i
				j_star = j
			end
		end
	end

	return (i_star, j_star)
end

iris = Iris()

# Normalize inputs
norm_const = maximum([maximum(iris.dataframe[:, i]) for i in axes(iris.features, 2)])
features_norm = iris.features ./ norm_const
setosa = filter(row -> occursin("setosa", row.class), iris.dataframe)[!, 1:4]
versi = filter(row -> occursin("versicolor", row.class), iris.dataframe)[!, 1:4]
virgi = filter(row -> occursin("virginica", row.class), iris.dataframe)[!, 1:4]

setosa_norm = setosa ./ norm_const
versi_norm = versi ./ norm_const
virgi_norm = virgi ./ norm_const

test = CSV.read("src/neural-networks/kohonen-iris/test_data.csv", DataFrame)
test_norm = test ./ norm_const

Wdf = CSV.read("src/neural-networks/kohonen-iris/model_weights.csv", DataFrame)

# FIXME: This script currently only works for a 50x50 model
W = Matrix(Wdf)
W = reshape([W[:, i] for i in 1:size(W, 2)], 50, 50)

setosa_x, versi_x, virgi_x, test_x = [], [], [], []
setosa_y, versi_y, virgi_y, test_y = [], [], [], []
for row in Vector.(eachrow(setosa_norm))
	coords = find_nearest(row, W)
	push!(setosa_x, coords[1])
	push!(setosa_y, coords[2])
end
for row in Vector.(eachrow(versi_norm))
	coords = find_nearest(row, W)
	push!(versi_x, coords[1])
	push!(versi_y, coords[2])
end
for row in Vector.(eachrow(virgi_norm))
	coords = find_nearest(row, W)
	push!(virgi_x, coords[1])
	push!(virgi_y, coords[2])
end
for row in Vector.(eachrow(test_norm))
	coords = find_nearest(row, W)
	push!(test_x, coords[1])
	push!(test_y, coords[2])
end

total_x = append!([], setosa_x, versi_x, virgi_x)
total_y = append!([], setosa_y, versi_y, virgi_y)
test_vals = collect(zip(test_x, test_y))

function plot_umatrix(W::Array)
	distances = zeros(size(W, 1), size(W, 2))
	for i in axes(W, 1)
		for j in axes(W, 2)
			for m in -1:1
				if i + m < 1 || i + m > size(W, 1)
					continue
				end
				for l in -1:1
					if j + l < 1 || j + l > size(W, 2)
						continue
					end
					distances[i, j] += norm(W[i, j] - W[i+m, j+l])
				end
			end
		end
	end

	return heatmap(
		distances,
		title="Overlay of data points on U-Matrix",
		size=(600, 550),
	)
end

function plot_analysis(W, xres=1000, yres=800, fontsize=8)
	weight_magnitudes = [norm(W[i, j]) for i in axes(W, 1), j in axes(W, 2)]
	sepal_area =
		[W[i, j][1] * W[i, j][2] * norm_const^2 for i in axes(W, 1), j in axes(W, 2)]
	petal_area =
		[W[i, j][3] * W[i, j][4] * norm_const^2 for i in axes(W, 1), j in axes(W, 2)]

	return plot(
		heatmap(weight_magnitudes),
		plot_umatrix(W),
		contourf(1:50, 1:50, petal_area),
		contourf(1:50, 1:50, sepal_area),
		layout=(2, 2),
		legend=false,
		size=(xres, yres),
		title=["Magnitude of weights" "U-Matrix" "Petal area" "Sepal area"],
		legendfontsize=fontsize,
		titlefontsize=round(Int, 1.5fontsize),
	)
end

function plot_overlay(type::Symbol, W, xres=1000, yres=800, markersize=3, fontsize=8)
	if type == :umatrix || type == :umat
		umat = plot_umatrix(W)
		umat = scatter!(
			[setosa_y, versi_y, virgi_y],
			[setosa_x, versi_x, virgi_x],
			label=["Setosa" "Versicolor" "Virginica"],
			title="Data points on closest neuron",
			marker=[:o :diamond :hexagon],
			markercolor=[:orange :white :green1],
			markersize=markersize,
			legend=:outerbottom,
			legend_column=-1,
			size=(xres, yres),
			legendfontsize=fontsize,
			titlefontsize=round(Int, 1.5fontsize),
		)

		# Test
		umat = scatter!(
			test_y,
			test_x,
			label="Test",
			marker=:rect,
			markercolor=:cyan,
			markersize=markersize,
		)

		return umat
	end

	if type == :all
		y = [setosa_y, versi_y, virgi_y]
		x = [setosa_x, versi_x, virgi_x]
		p = plot_analysis(W)
		p = scatter!(
			[y, y, y, y],
			[x, x, x, x],
			label=["Setosa" "Versicolor" "Virginica"],
			marker=[:o :diamond :hexagon],
			markercolor=[:orange :white :green1],
			markersize=markersize,
			legend=:outerbottom,
			legend_column=-1,
			size=(xres, yres),
			legendfontsize=fontsize,
			titlefontsize=round(Int, 1.5fontsize),
			colorbar_tickfontsize=fontsize,
		)

		# Test
		p = scatter!(
			[test_y, test_y, test_y, test_y],
			[test_x, test_x, test_x, test_x],
			label="Test",
			marker=:rect,
			markercolor=:cyan,
			markersize=markersize,
		)

		return p
	end
end

# TODO: Fix this up and move it to a different file
function train_predictor(; quiet=false)
	n = 2
	ϵ = 0.2
	E_stop = 10^-3
	max_iter = 100_000

	examples = hcat(total_x, total_y)
	norm_const = maximum(examples)
	examples = 2 * (examples ./ norm_const) .- 1
	#              Setosa = 1   Versicolor = 0 Virginica = -1
	results = vcat(fill(1, 50), fill(0, 50), fill(-1, 50))

	yarr = rand((-1, 1)) * rand(length(results))
	w = rand((-1, 1)) * rand(n)

	cycle = 0
	while cycle < max_iter
		ν = rand(1:length(results))
		x_ν = examples[ν, :]
		y = tanh(x_ν ⋅ w)
		yarr[ν] = y

		δ = results[ν] - y
		w += ϵ * δ * (1 - tanh(y)^2) * x_ν

		E = 0.5 * sum([(ξ_ν - y_ν)^2 for (ξ_ν, y_ν) in zip(results, yarr)])
		if cycle % 1000 == 0 && !quiet
			println("Cycle $cycle. Error: $E")
		end
		if E < E_stop
			break
		end
		cycle += 1
	end

	return w
end

classify(x::Vector{Tuple{Any, Any}}, w::Vector) = [classify(x_i, w) for x_i in x]
classify(x::Vector{Vector}, w::Vector) = [classify(x_i, w) for x_i in x]
classify(x::Tuple, w::Vector) = classify(collect(x), w)
function classify(x::Vector, w::Vector)
	y = tanh(x ⋅ w)

	d1 = abs(1 - y)
	d2 = abs(y)
	d3 = abs(-1 - y)
	if d1 < d2 && d1 < d3
		return "Iris Setosa"
	elseif d2 < d1 && d2 < d3
		return "Iris Versicolor"
	elseif d3 < d1 && d3 < d2
		return "Iris Virginica"
	end
end