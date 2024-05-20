using Amaru

fc = -30e3
ft = 3e3
fb = -33e3

fc_fun = PathFunction(:M, 0.0, 0.4*fc, :L, 0.003, 1.0*fc, :L, 0.005, 0.2*fc)
ft_fun = PathFunction(:M, 0.0, 1.0*ft, :L, 0.0002, 0)

fic_fun = PathFunction(:M, 0.0, 0.66*fc, :L, 0.002, 1.5*fc)


# Plotting

chart = Chart(
)

X = range(0, 0.006, 200)
Y = fic_fun.(X)

series = [ 
    LineSeries(X, Y),
]

addseries!(chart, series)
save(chart, "fun.pdf")