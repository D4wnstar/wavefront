using LinearAlgebra, DataFrames, CSV, MLDatasets, Plots, Random

iris = Iris()

# Plotting data
setosa = iris.dataframe[iris.dataframe.class.=="Iris-setosa", :]
versi = iris.dataframe[iris.dataframe.class.=="Iris-versicolor", :]
virgi = iris.dataframe[iris.dataframe.class.=="Iris-virginica", :]

# Dataset plots
function make_plots(type::Symbol)
	if type == :petals
		p = scatter(
			[
				setosa.petallength,
				versi.petallength,
				virgi.petallength,
			],
			[
				setosa.petalwidth,
				versi.petalwidth,
				virgi.petalwidth,
			],
			label=["Setosa" "Versicolor" "Virginica"],
			xlims=(0.9, 7.1),
			ylims=(0, 2.6),
			xlabel="Petal length [cm]",
			ylabel="Petal width [cm]",
			marker=[:o :diamond :hexagon],
		)
		display(p)
	end

	if type == :sepals
		p = scatter(
			[
				setosa.sepallength,
				versi.sepallength,
				virgi.sepallength,
			],
			[
				setosa.sepalwidth,
				versi.sepalwidth,
				virgi.sepalwidth,
			],
			label=["Setosa" "Versicolor" "Virginica"],
			xlims=(4.0, 8.0),
			ylims=(1.9, 4.6),
			xlabel="Sepal length [cm]",
			ylabel="Sepal width [cm]",
			marker=[:o :diamond :hexagon],
		)
		display(p)
	end
end

# Heatmap for weights
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
					distances[i, j] += norm(W[i, j, :] - W[i+m, j+l, :])
				end
			end
		end
	end

	return heatmap(distances)
end

# Training
norm_const = maximum([maximum(iris.dataframe[:, i]) for i in axes(iris.features, 2)])
features = iris.features ./= norm_const
petals = features[!, r"petal"]
sepals = features[!, r"sepal"]

dsq(i, i_star, j, j_star) = (i - i_star)^2 + (j - j_star)^2
θ(i, i_star, j, j_star, σ) = exp(-(dsq(i, i_star, j, j_star) / 2σ))

function find_nearest(x::Vector, W::Array)
	i_star = 1
	j_star = 1
	shortest_dist = Inf
	for i in axes(W, 1)
		for j in axes(W, 2)
			curr_dist = norm(x - W[i, j, :])
			if curr_dist < shortest_dist
				shortest_dist = curr_dist
				i_star = i
				j_star = j
			end
		end
	end

	return (i_star, j_star)
end

function train_model(
	data_name::Symbol,
	h,
	k;
	ϵ_stop=0.001,
	max_cycle=1_000_000,
	filename="model_weights.csv",
)
	if data_name == :petals
		n_inputs = 2
		data = petals
	elseif data_name == :sepals
		n_inputs = 2
		data = sepals
	elseif data_name == :all
		n_inputs = 4
		data = DataFrame([sepals petals])
	else
		error("Invalid data type. Must be :petals, :sepals or :all")
	end

	Random.seed!(123)
	W = rand(h, k, n_inputs)
	measurements = size(data, 1)
	ϵ = 0.2
	σ = 16
	cycle = 1
	training = Animation()

	while cycle <= max_cycle && ϵ > ϵ_stop
		ν = data[rand(1:measurements), :]
		x = Vector(ν)
		i_star, j_star = find_nearest(x, W)

		for i in axes(W, 1)
			for j in axes(W, 2)
				W[i, j, :] .+= ϵ * θ(i, i_star, j, j_star, σ) * (x - W[i, j, :])
			end
		end

		if cycle % 100 == 0
			ϵ = 0.995ϵ
			σ = 0.98σ
		end

		if (cycle == 1 || cycle % 1000 == 0)
			println("Cycle $cycle. ϵ: $ϵ. σ: $σ")

			if n_inputs == 2
				lengths, widths = [], []
				for i in axes(W, 1)
					for j in axes(W, 2)
						# Add round(..., digits=1) on the next two values to make the neurons snap to a grid
						l = W[i, j, 1] * norm_const
						w = W[i, j, 2] * norm_const
						push!(lengths, l)
						push!(widths, w)
					end
				end

				if data_name == :petals
					object = "Petal"
					xlims = (0.9, 7.1)
					ylims = (0, 2.6)
				else
					object = "Sepal"
					xlims = (4.0, 8.0)
					ylims = (1.9, 4.6)
				end

				frame(
					training,
					scatter(
						lengths,
						widths,
						xlabel="$object length [cm]",
						ylabel="$object width [cm]",
						xlims=xlims,
						ylims=ylims,
						legend=false,
					),
				)
			end
		end

		cycle += 1
	end

	if n_inputs == 2
		gif(training, "src/neural-networks/kohonen-iris/training_montage.gif", fps=10)
	end

	weight_magnitudes = [norm(W[i, j, :]) for i in axes(W, 1), j in axes(W, 2)]
	p = plot(
		heatmap(weight_magnitudes),
		plot_umatrix(W),
		layout=(1, 2),
		legend=false,
		size=(1000, 400),
		title=["Magnitude of weights" "U-Matrix"],
	)
	display(p)

	Wdf = DataFrame()
	for i in axes(W, 1)
		for j in axes(W, 2)
			Wdf[!, "$i-$j"] = W[i, j, :]
		end
	end

	CSV.write("src/neural-networks/kohonen-iris/$filename", Wdf)
end
