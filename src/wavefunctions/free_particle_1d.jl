using LinearAlgebra
using PlotThemes
using Plots
using LaTeXStrings
using QuadGK

theme(:dracula)

# ħ = 1.0545e-34#Js
ħ = 2

ψ_0(y, p0) = (2 / pi)^(1 / 4) * exp(p0^2 / 4ħ^2) * exp(-y^2 + im * p0 * y / ħ)

propagator(x, y, t, m) = sqrt(m / (2im * t * π * ħ)) * exp(im * m * (x - y)^2 / (2ħ * t))

function ψ_t(x, t, m, p0)
	quadgk(y -> propagator(x, y, t, m) * ψ_0(y, p0), -Inf, +Inf, rtol=0.0001)[1]
end

function make_alphas(len::Integer; decay=0.1)
	alphas = ones(len)
	for i in eachindex(alphas)
		alphas[i] = max(alphas[i] - decay * (len - i), 0)
	end
	return alphas
end

x = -2:0.1:5
m = 2
p_0 = 1

y_vals = []
alphas = []

@time @gif for (i, t) in enumerate(0:0.05:4)
	println("Calculating time $t...")
	if i == 1
		init_state = norm.(ψ_0.(x, p_0)) .^ 2
		push!(y_vals, init_state)
		alphas = [1]'
	else
		new_state = @. norm(ψ_t(x, t, m, p_0))^2
		push!(y_vals, new_state)
		alphas = make_alphas(i, decay=0.2)'
	end

	plot(
		x,
		y_vals,
		legend=false,
		alpha=alphas,
		xlims=(-2, 5),
		ylims=(0, 1),
		xlabel=L"x",
		ylabel=L"|\psi(x,t)|^2",
		color=:cyan,
		title=L"Free particle position distribution at time $t=%$t$",
	)
end fps = 12
