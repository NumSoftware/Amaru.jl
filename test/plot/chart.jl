using Amaru

X = collect(0:0.4:8)
Y1 = sin.(X)
Y2 = cos.(X)
Y3 = sin.(X).*cos.(X)

chart = Chart(
    x_label = L"$x$ coordinate",
    y_label = L"$y$ coordinate",
    legend = :bottomright
)

series = [
    LineSeries(X, Y1, marker=:circle, lw=0.5, label=L"\sin(x)")
    LineSeries(X, Y2, marker=:utriangle, lc=:royalblue, label=L"\cos(x)")
    LineSeries(X, Y3, marker=:square, lc=:green, ls=:dash, label=L"\sin(x)\, \cos(x)")
]

addseries!(chart, series)
addlegend!(chart, Legend(ncols=3, loc=:outerbottom))
addannotation!(chart, Annotation(L"Amaru $\frac{x}{y}$", 0.5, 0.9, target=(1.5,1), textalignment=:top))
addannotation!(chart, Annotation(L"Amaru $\frac{x}{y}$", 0.7, 0.2, target=(4,0), textalignment=:top))
# addannotation!(chart, Annotation(L"Amaru $\frac{x}{y}$", 0.2, 0.7, target=(3.5,0), textalignment=:top))
# addannotation!(chart, Annotation(L"Amaru $\frac{x}{y}$", 0.55, 0.9, target=(3.5,0), textalignment=:top))

save(chart, "chart.pdf")