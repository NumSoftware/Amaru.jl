using Amaru

X = collect(0:0.4:8)
Y1 = sin.(X)
Y2 = cos.(X)

chart = Chart(
    xlabel = L"$x$ coordinate",
    ylabel = L"$y$ coordinate",
    legend = :bottomright
)

p1 = LinePlot(X, Y1, marker=:circle, lw=0.5, label="sin")
p2 = LinePlot(X, Y2, marker=:utriangle, lc=:royalblue, ls=:dash, label="cos")

addplot!(chart, p1, p2)

save(chart, "chart.pdf")