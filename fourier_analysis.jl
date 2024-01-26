using Plots
using DynamicQuantities
using DynamicQuantities.Units: mm, cm, nm, m, μm
using FFTW

x = 0μm:0.1μm:4π*μm
wavefunc(x) = sin(x)
y = wavefunc.(x)
transform = fft(y)

plot(x, [y transform], layout=(2, 1))