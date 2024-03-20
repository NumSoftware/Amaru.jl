using Amaru

X = collect(0:0.4:8)
Y1 = sin.(X)
Y2 = cos.(X)
Y3 = sin.(X).*cos.(X)

chart = Chart(
    xlabel = L"$x$ coordinate",
    ylabel = L"$y$ coordinate",
    legend = :bottomright
)

series = [
    LineSeries(X, Y1, marker=:circle, lw=0.5, label=L"\sin(x)")
    LineSeries(X, Y2, marker=:utriangle, lc=:royalblue, label=L"\cos(x)")
    LineSeries(X, Y3, marker=:square, lc=:green, ls=:dash, label=L"\sin(x)\, \cos(x)")
]

addseries!(chart, series)
addlegend!(chart, Legend(ncols=3, loc=:outertop))

save(chart, "chart.pdf")