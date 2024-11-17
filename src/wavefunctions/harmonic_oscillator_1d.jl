using LinearAlgebra
using PlotThemes
using Plots
using LaTeXStrings

theme(:dracula)

# ħ = 1.0545e-34#Js
ħ = 2

# Excited states have been derived from Hermite functions in natural units
# then SI units have been put back in and renormalized
ψ_0(x, m, ω) = (m * ω / (π * ħ))^(1 / 4) * exp(-(m * ω / 2ħ) * x^2)
ψ_1(x, m, ω) = (m * ω / (π * ħ))^(1 / 4) * sqrt(2m * ω / ħ) * x * exp(-(m * ω / 2ħ) * x^2)

x = -3:0.01:3
m = 2
ω = 1
p_0 = 1

y0 = @. norm(ψ_0(x, m, ω))^2
y1 = @. norm(ψ_1(x, m, ω))^2

plot(
	x,
	[y0 y1],
	label=["Ground state" "First excited state"],
	ylims=(0, 1),
	xlabel=L"x",
	ylabel=L"|\psi|^2",
	title="Quantum harmonic oscillator position distribution for eigenstates",
	size=(800, 600),
)
