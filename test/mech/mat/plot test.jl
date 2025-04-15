using Amaru

fc = -30e3
ft = 3e3
fb = -33e3

fc_fun = PathFunction(:M, 0.0, 0.4*fc, :L, 0.003, 1.0*fc, :L, 0.005, 0.2*fc)
fc_fun = PathFunction(:M, 0.0, 0.4*fc, :Q, 0.003, fc, 0.02, fc)
fc_fun = PathFunction(:M, 0.0, 0.4*fc, :Q, 0.003, 0.4*fc, 0.02, fc)

ft_fun = PathFunction(:M, 0.0, ft, :L, 0.0002, 0)
ft_fun = PathFunction(:M, 0.0, ft, :C, 0.0005, 0.1*ft, 0.002, 0, 0.003, 0 )

fic_fun = PathFunction(:M, 0.0, 0.66*fc, :L, 0.002, 1.5*fc)

# Plotting

chart = Chart()

X = range(0, 0.02, 200)
# X = range(0, 0.003, 200)
Y = fc_fun.(X)



series = [ 
    LineSeries(X, -Y),
]

addseries!(chart, series)
save(chart, "fun.pdf")
