using Plots
using FFTW

# As the diffraction pattern of a slit is effectively the (absolute value of the) Fourier
# transform of the electric field incident on the slit over the surface of the slit,
# taking the transform of the diffraction pattern effectively undoes the diffraction
# and gives the shape of the slit
"Take a 2D diffraction pattern and render the shape of the slit from which it came."
function render_slit_shape(
	x::AbstractVector,
	y::AbstractVector,
	diff_pattern::AbstractMatrix;
	zoom::Int=10,
)
	transform = fft(diff_pattern) |> fftshift
	abs_transform = abs.(transform)
	zoom_limit_x, zoom_limit_y = size(abs_transform) .÷ zoom
	x_half = length(x) ÷ 2
	y_half = length(y) ÷ 2
	heatmap(
		x,
		y,
		[abs_transform, abs_transform],
		title=["Slit shape from diffraction pattern" "Zoomed in"],
		aspectratio=1,
		xlims=[:auto (x[x_half-zoom_limit_x], x[x_half+zoom_limit_x])],
		ylims=[:auto (y[y_half-zoom_limit_y], y[y_half+zoom_limit_y])],
		layout=(1, 2),
		size=(800, 400),
	)
end
