export newchart, savechart, showchart, cplot

const CHART_ARGS = Dict{String, Any}()

"""
    newchart(kwargs...)

Configures the ploting area for a `PyPlot` line chart.

# Keyword arguments

`xlabel         = "\$x\$"` : x label

`ylabel         = "\$y\$"` : y label

`lw             = 0.7` : line width

`xbins          = 6` : number of bins in x

`ybins          = 6` : number of bins in y

`grid           = false` : grid

`figsize        = (3,2)` : size of plot

`fontsize       = 7` : font size

`legendloc      = "best"` : legend location in plot

`legendfontsize = 0` : defaults to fontsize

`legendexpand   = false` : expanded version

`labelspacing   = 0.5` : spacing between legend labels

`ticksinside    = true` : put ticks inside plot

# Example

The following example shows the function usage:

```
using Amaru, PyPlot

newchart()
X = 1:10
Y = rand(10)
Z = rand(10)

plot(X, Y, "k", label="label1")
plot(X, Z, "b", label="label2")

savechart("filename.pdf")
```
"""
function newchart(; 
                  xlabel         = "\$x\$",
                  ylabel         = "\$y\$",
                  lw             = 0.7,
                  xbins          = 6,
                  ybins          = 6,
                  grid           = false,
                  figsize        = (3,2),
                  fontsize       = 7,
                  legendloc      = "best",
                  legendfontsize = 0,
                  legendexpand   = false,
                  labelspacing   = 0.5,
                  ticksinside    = true
              )

    CHART_ARGS["legendloc"]    = legendloc
    CHART_ARGS["legendexpand"] = legendexpand
    CHART_ARGS["labelspacing"] = labelspacing

    @eval import PyPlot:plt,gca

    # Configure plot
    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", scale_dashes=true)

    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("lines", lw=0.7)
    plt.rc("lines", markersize=2)
    plt.rc("axes" , linewidth=0.5)
    plt.rc("figure", figsize=(3, 2))
    legendfontsize==0 && (legendfontsize=fontsize)
    plt.rc("legend", fontsize=legendfontsize)

    grid && plt.grid( color="lightgrey", ls="dotted", lw=0.3)
    ax = plt.axes()
    plt.locator_params(axis="x", nbins=xbins)
    plt.locator_params(axis="y", nbins=ybins)

    # Set limits
    !isnan(xmax) && plt.xlim(right=xmax)
    !isnan(ymax) && plt.ylim(top=ymax)

    # Print axes labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Tick marks direction
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis="x", direction="in")

    return nothing
end


"""
    savechart(filename; kwargs...)

Saves the current `PyPlot` line chart.

# Arguments

`filename` : output file

# Keyword arguments

`showlegend = true` : controls if legend is shown

`format     = "pdf"` : output file format

"""
function savechart(filename; showlegend=true, format="pdf")
    legendloc    = CHART_ARGS["legendloc"]
    legendexpand = CHART_ARGS["legendexpand"]
    labelspacing = CHART_ARGS["labelspacing"]
    plt.ioff()
    ax = gca()
    if showlegend
        mode = legendexpand ? "expand" : nothing
        ncol = length(ax.lines)

        if legendloc=="top"
            leg = plt.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), edgecolor="k", ncol=ncol, mode=mode)
        elseif legendloc=="right"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1), edgecolor="k")
        elseif legendloc=="bottom"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.02, 1.04, -0.2), edgecolor="k", ncol=ncol, mode=mode)
        else
            leg = plt.legend(loc=legendloc, edgecolor="k", labelspacing=labelspacing)
        end

        frame = leg.get_frame()
        frame.set_linewidth(0.5)
    end
    plt.savefig(filename, bbox_inches="tight", pad_inches=0.01, format=format)
end


function showchart(showlegend=true, legendloc="best", labelspacing=0.5, legendexpand=false)
    plt.ion()
    plt.legend(loc=legendloc, edgecolor="k", labelspacing=labelspacing)
    show()
end


function cplot(X, Y, filename=""; xlabel="\$x\$", ylabel="\$y\$", lw=0.7, ls="-", ms=2, marker=nothing, color="", legend=[], legendloc="best",
               xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
               xmin=NaN, xmax=NaN,
               ymin=NaN, ymax=NaN,
               fontsize=7, legendfontsize=0, labelspacing=0.5)

    @eval import PyPlot:plt, matplotlib, figure

    # Fix labels
    texlabels = Dict(
                  "sxx"=>raw"$\sigma_{xx}$",
                  "syy"=>raw"$\sigma_{yy}$",
                  "szz"=>raw"$\sigma_{zz}$",
                  "sxy"=>raw"$\sigma_{xy}$",
                  "syz"=>raw"$\sigma_{yz}$",
                  "sxz"=>raw"$\sigma_{xz}$",
                  "exx"=>raw"$\varepsilon_{xx}$",
                  "eyy"=>raw"$\varepsilon_{yy}$",
                  "ezz"=>raw"$\varepsilon_{zz}$",
                  "exy"=>raw"$\varepsilon_{xy}$",
                  "eyz"=>raw"$\varepsilon_{yz}$",
                  "exz"=>raw"$\varepsilon_{xz}$",
                  "sigma" =>raw"$\sigma$",
                  "tau"   =>raw"$\tau$",
                  "lambda"=>raw"$\lambda$",
                 )

    labels = [ xlabel, ylabel ]
    for (i,label) in enumerate(labels)
        if haskey(texlabels, label)
            label = get(texlabels, label, label)
        else
            !occursin(r"^\$.+\$$", label) && ( label = replace(label, r"(\w+)_(\w+)" =>  s"$\1_{\2}$") )
            if occursin(r"^\$.+\$$", label)
                label = replace(label, r"\b(sigma|tau|lambda)(\b|_)" => s"\\\1\2")
            end
        end
        labels[i] = label
    end
    xlabel, ylabel = labels

    # Configure plot
    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", scale_dashes=true)

    if filename!=""
        plt.ioff()
        plt.rc("xtick", labelsize=7)
        plt.rc("ytick", labelsize=7)
        plt.rc("lines", lw=0.7)
        plt.rc("lines", markersize=2)
        plt.rc("axes" , linewidth=0.5)
        plt.rc("figure", figsize=(3, 2))
        legendfontsize==0 && (legendfontsize=fontsize)
        plt.rc("legend", fontsize=legendfontsize)
    end

    # Convert amtrices to array of vectors
    size(X,2)>1 && ( X = [ X[:,j] for j=1:size(X,2) ] )
    size(Y,2)>1 && ( Y = [ Y[:,j] for j=1:size(Y,2) ] )

    # Convert vectors to array of one vector
    eltype(X)<:Real && (X = [X] )
    eltype(Y)<:Real && (Y = [Y] )
    mx = length(X)
    my = length(Y)
    mx == 1 && ( X = repeat(X, my) )
    my == 1 && ( Y = repeat(Y, mx) )
    @assert length(X)==length(Y) "cplot: X and Y should have the same size"
    m = max(mx, my)

    # Fix markers
    if typeof(marker)==String || marker==nothing
        marker = repeat([marker], m)
    end
    @assert length(marker)==m

    # Fix markers size
    if typeof(ms) <: Real
        ms = repeat([ms], m)
    end
    @assert length(ms)==m

    # Fix line widths - lw
    if typeof(lw) <: Real
        lw = repeat([lw], m)
    end
    @assert length(lw)==m

    # Fix line styles - ls
    if typeof(ls)==String
        ls = repeat([ls], m)
    end
    @assert length(ls)==m

    # Fix line colors
    if typeof(color)==String
        if color==""
            color = nothing
        end
        color = repeat([color], m)
    end
    @assert length(color)==m

    # Fix legend labels
    haslegend = true
    if length(legend)==0
        haslegend = false
        legend = repeat([""],m)
    end
    @assert length(legend)==m


    # Plot curves
    for i=1:m
        plt.plot(X[i], Y[i], marker=marker[i], ms=ms[i], color=color[i], lw=lw[i], ls=ls[i], label=legend[i])
    end

    grid && plt.grid( color="lightgrey", ls="dotted", lw=0.3)
    ax = plt.axes()
    plt.locator_params(axis="x", nbins=xbins)
    plt.locator_params(axis="y", nbins=ybins)

    # Set limits
    !isnan(xmax) && plt.xlim(right=xmax)
    !isnan(ymax) && plt.ylim(top=ymax)

    # Print axes labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # plot legend
    if haslegend
        mode = legendexpand ? "expand" : nothing
        ncol = ncol==0 ? m : ncol

        if legendloc=="top"
            leg = plt.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), edgecolor="k", ncol=ncol, mode=mode)
        elseif legendloc=="right"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1), edgecolor="k")
        elseif legendloc=="bottom"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.02, 1.04, -0.2), edgecolor="k", ncol=ncol, mode=mode)
        else
            leg = plt.legend(loc=legendloc, edgecolor="k", labelspacing=labelspacing)
        end

        frame = leg.get_frame()
        frame.set_linewidth(0.5)
    end

    # show or save plot
    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.01, format="pdf")
    end

    plt.close("all")
end


#TODO
using LaTeXStrings
# cplot([
#    (x=X, y=Y, lw=1, ls="-", ms=2, marker="s", color="", label="")
#    (x=X, y=Y, lw=1, ls="-", ms=2, marker="s", color="", label="")
#    (x=X, y=Y, lw=1, ls="-", ms=2, marker="s", color="", label="")
#    (x=X, y=Y, lw=1, ls="-", ms=2, marker="s", color="", label="")
#    ],
#    out.pdf,
#    xlabel=L"x", ylabel=L"y", legendloc="best",
#    xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
#    xmin=NaN, xmax=NaN,
#    ymin=NaN, ymax=NaN,
#    fontsize=7, legendfontsize=0, labelspacing=0.5)
#
#=
function cplot2(
                   data::Array{<:NamedTuple}, 
                   filename::String="";
                   xlabel=L"x", 
                   ylabel=L"y", 
                   legendloc="best",
                   xbins=6, 
                   ybins=6, 
                   grid=false, 
                   figsize=(3,2), 
                   legendexpand=false, 
                   ncol=0,
                   xmin=NaN, 
                   xmax=NaN, 
                   ymin=NaN, 
                   ymax=NaN,
                   fontsize=7, 
                   legendfontsize=0, 
                   labelspacing=0.5, 
                   backend=:pyplot
                  )

    kwdargs = (xlabel=xlabel, ylabel=ylabel, legendloc=legend, xbins=xbins, ybins=ybins, 
               grid=grid, figsize=figsize, legendexpand=legendexpand, ncol=ncol,
               xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fontsize=fontsize, 
               legendfontsize=legendfontsize, labelspacing=labelspacing, backend=backend)


    #backend = get(args, :backend, :pyplot)
    if backend==:pyplot
        cplot_pyplot(data, filename; kwdargs...)
    elseif backend==:pgf
        cplot_pgf(data, filename; kwdargs...)
    else
        error("wrong backend")
    end
end

function cplot_pgf(
                  )
end
=#

function cplot_pyplot(data::Array{<:NamedTuple}, filename::String=""; 
                xlabel=L"x", ylabel=L"y", legendloc="best",
                xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
                xmin=NaN, xmax=NaN, ymin=NaN, ymax=NaN,
                fontsize=7, legendfontsize=0, labelspacing=0.5, backend=:pyplot)

    @eval import PyPlot:plt, matplotlib, figure

    # Fix labels
    texlabels = Dict(
                  "sxx"=>raw"$\sigma_{xx}$",
                  "syy"=>raw"$\sigma_{yy}$",
                  "szz"=>raw"$\sigma_{zz}$",
                  "sxy"=>raw"$\sigma_{xy}$",
                  "syz"=>raw"$\sigma_{yz}$",
                  "sxz"=>raw"$\sigma_{xz}$",
                  "exx"=>raw"$\varepsilon_{xx}$",
                  "eyy"=>raw"$\varepsilon_{yy}$",
                  "ezz"=>raw"$\varepsilon_{zz}$",
                  "exy"=>raw"$\varepsilon_{xy}$",
                  "eyz"=>raw"$\varepsilon_{yz}$",
                  "exz"=>raw"$\varepsilon_{xz}$",
                  "sigma" =>raw"$\sigma$",
                  "tau"   =>raw"$\tau$",
                  "lambda"=>raw"$\lambda$",
                 )

    if backend==:pyplot
    elseif backend==:pgfplots
    end

    # Configure plot
    plt.close("all")

    #=
    import matplotlib.font_manager as font_manager

    font_dirs = ['/my/custom/font/dir', ]
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    font_list = font_manager.createFontList(font_files)
    font_manager.fontManager.ttflist.extend(font_list)

    mpl.rcParams['font.family'] = 'My Custom Font'
    =#


    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", scale_dashes=true)

    if filename!=""
        plt.ioff()
        plt.rc("xtick", labelsize=7)
        plt.rc("ytick", labelsize=7)
        plt.rc("lines", lw=0.7)
        plt.rc("lines", markersize=2)
        plt.rc("axes" , linewidth=0.5)
        plt.rc("figure", figsize=(3, 2))
        legendfontsize==0 && (legendfontsize=fontsize)
        plt.rc("legend", fontsize=legendfontsize)
        #plt.rc("xtick.major.size") = 20
        #plt.rc("xtick.major.width") = 4
        #plt.rc("xtick.minor.size") = 10
        #plt.rc("xtick.minor.width") = 2
    end

    # Plot curves
    for (i, line) in enumerate(data)
        X = line.X
        Y = line.Y
        lw = get(line, :lw, 0.7)
        ms = get(line, :ms, 2)
        marker = get(line, :marker, nothing)
        color = get(line, :color, "")
        label = get(line, :label, "")
        plt.plot(X, Y, marker=marker, ms=ms, color=color, lw=lw, ls=ls, label=label)
    end

    grid && plt.grid( color="lightgrey", ls="dotted", lw=0.3)
    ax = plt.axes()
    plt.locator_params(axis="x", nbins=xbins)
    plt.locator_params(axis="y", nbins=ybins)

    # Set limits
    !isnan(xmax) && plt.xlim(right=xmax)
    !isnan(ymax) && plt.ylim(top=ymax)

    # Print axes labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # plot legend
    if haslegend
        mode = legendexpand ? "expand" : nothing
        ncol = ncol==0 ? m : ncol

        if legendloc=="top"
            leg = plt.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), edgecolor="k", ncol=ncol, mode=mode)
        elseif legendloc=="right"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1), edgecolor="k")
        elseif legendloc=="bottom"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.02, 1.04, -0.2), edgecolor="k", ncol=ncol, mode=mode)
        else
            leg = plt.legend(loc=legendloc, edgecolor="k", labelspacing=labelspacing)
        end

        frame = leg.get_frame()
        frame.set_linewidth(0.5)
    end

    # Tick marks direction
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(axis="y", direction="in")
    #ax.xaxis.set_tick_params(width=5)
    #ax.yaxis.set_tick_params(width=5)

    # show or save plot
    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.01, format="pdf")
    end

    plt.close("all")

end
