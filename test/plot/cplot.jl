# Loading the finite element package Amaru
using Amaru, LaTeXStrings
using Test


printstyled("\nPlotting\n", color=:blue, bold=true)

# 2D
X = range(0, pi/2, 10)
Y1 = sin.(X)
Y2 = cos.(X)
Y3 = sin.(X).^2

cplot( 
    (x=X, y=Y1, marker="o", label=L"\sin(x)"),
    (x=X, y=Y2, ls="--", label=L"\cos(x)"),
    (x=X, y=Y3, ls="--", tagpos=0.5, tag=L"\sin^2(x)"),
    "out.pdf",
)
println( @test isfile("out.pdf") )
rm("out.pdf")

cplot( 
    (x=X, y=Y1, label=L"\sin(x)"),
    (x=X, y=Y2, label=L"\cos(x)"),
    "out.pdf",
    xlabel    = L"x",
    ylabel    = L"y",
    xbins     = 7,
    ybins     = 5,
    xmult     = 10.0,
    ymult     = 2.0,
    showgrid  = false,
    xlim      = (0,10),
    legendloc = 6,
    figsize   = (3,2),
    quiet     = false
)
println( @test isfile("out.pdf") )
rm("out.pdf")
