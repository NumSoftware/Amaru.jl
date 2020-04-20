using LaTeXStrings
export newchart, savechart, showchart, cplot, cplot0

#const CHART_ARGS = Dict{String, Any}()

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

    #CHART_ARGS["legendloc"]    = legendloc
    #CHART_ARGS["legendexpand"] = legendexpand
    #CHART_ARGS["labelspacing"] = labelspacing

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


function cplot0(X, Y, filename=""; xlabel="\$x\$", ylabel="\$y\$", lw=0.7, ls="-", ms=2, marker=nothing, color="", legend=[], legendloc="best",
               xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
               xlim = nothing,
               ylim = nothing,
               #xmin=NaN, xmax=NaN,
               #ymin=NaN, ymax=NaN,
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



"""
    cplot(data, filename, kwargs...)

Plots a pyplot line chart.

# Arguments

`data` : An array of named tuples with information for each line

`filename` = "" : The chart filename

# Keyword arguments

`xlabel         = "\$x\$"` : x label

`ylabel         = "\$y\$"` : y label

`legendloc      = "best"` : legend location in plot

`xbins          = 6` : number of bins in x

`ybins          = 6` : number of bins in y

`grid           = false` : grid

`figsize        = (3,2)` : size of plot

`legendexpand   = false` : expanded version

`ncol           = 0` : number of columns in horizontal legend

`xlim           = nothing` : x axis limits

`ylim           = nothing` : y axis limits

`xscale         = "linear"` : x axis scale

`yscale         = "linear"` : y axis scale

`ticksinside    = true` : put ticks inside plot

`fontsize       = 7` : font size

`legendfontsize = 0` : defaults to fontsize

`labelspacing   = 0.5` : spacing between legend labels

# Example

The following example shows the function usage:

```
using Amaru
X1 = collect(1:20)
X2 = collect(1:20)
Y1 = rand(20)
Y2 = rand(20)

cplot([
   (x=X1, y=Y1, lw=1, ls="-", ms=2, marker="s", color="", label="")
   (x=X2, y=Y2, lw=1, ls="-", ms=2, marker="s", color="", label="")
   ],
   "out.pdf",
   xlabel=raw"\$x\$", ylabel="\$y\$", legendloc="best",
   xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
   xlim=(0,1), xscale="linear", yscale="log",
   fontsize=7, legendfontsize=0, labelspacing=0.5)
```
"""
function cplot(data::Array{<:NamedTuple}, 
               filename::String = "";
               xlabel           = L"x",
               ylabel           = L"y",
               legendloc        = "best",
               xbins            = 6,
               ybins            = 6,
               grid             = false,
               figsize          = (3,2),
               legendexpand     = false,
               ncol             = 0,
               xlim             = nothing,
               ylim             = nothing,
               xscale           = "linear",
               yscale           = "linear",
               ticksinside      = true,
               fontsize         = 7,
               legendfontsize   = 0,
               labelspacing     = 0.5)

    if filename==""
        printstyled("cplot: generating plot\n", color=:cyan )
    else
        printstyled("cplot: generating plot to file ", repr(filename), "\n", color=:cyan )
    end
    printstyled("Available options: xlabel, ylabel, legendloc, xbins, ybins, grid, figsize, legendexpand, ncol,"*
            " xlim, ylim, xscale, yscale, fontsize, ticksinside, legendfontsize, labelspacing\n", color=:light_black)
    printstyled("Options per curve: x, y, color, ls, lw, marker, ms, label\n", color=:light_black)


    @eval import PyPlot:plt, matplotlib, figure

    line_styles = ("-", "--", "-.", ":", "", " ", "None", nothing)
    markers     = (".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "*", "h", "H", "+", "x", "D", "d", "|", "_", "P", "X", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, "None", nothing, " ", "")
    colors      = ("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
                   "c", "b", "w", "g", "y", "k", "r", "m",
                   "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
                   "indigo", "gold", "hotpink", "firebrick", "indianred", "yellow", "mistyrose", "darkolivegreen", "olive", "darkseagreen", "pink", "tomato", "lightcoral", "orangered", "navajowhite", "lime", "palegreen", "darkslategrey", "greenyellow", "burlywood", "seashell", "mediumspringgreen", "fuchsia", "papayawhip", "blanchedalmond", "chartreuse", "dimgray", "black", "peachpuff", "springgreen", "aquamarine", "white", "orange", "lightsalmon", "darkslategray", "brown", "ivory", "dodgerblue", "per", "lawngreen", "chocolate", "crimson", "forestgreen", "darkgrey", "lightseagreen", "cyan", "mintcream", "silver", "antiquewhite", "mediumorchid", "skyblue", "gray", "darkturquoise", "goldenrod", "darkgreen", "floralwhite", "darkviolet", "darkgray", "moccasin", "saddlebrown", "grey", "darkslateblue", "lightskyblue", "lightpink", "mediumvioletred", "slategrey", "red", "deeppink", "limegreen", "darkmagenta", "palegoldenrod", "plum", "turquoise", "lightgrey", "lightgoldenrodyellow", "darkgoldenrod", "lavender", "maroon", "yellowgreen", "sandybrown", "thistle", "violet", "navy", "magenta", "dimgrey", "tan", "rosybrown", "olivedrab", "blue", "lightblue", "ghostwhite", "honeydew", "cornflowerblue", "slateblue", "linen", "darkblue", "powderblue", "seagreen", "darkkhaki", "snow", "sienna", "mediumblue", "royalblue", "lightcyan", "green", "mediumpurple", "midnightblue", "cornsilk", "paleturquoise", "bisque", "slategray", "darkcyan", "khaki", "wheat", "teal", "darkorchid", "salmon", "deepskyblue", "rebeccapurple", "darkred", "steelblue", "palevioletred", "lightslategray", "aliceblue", "lightslategrey", "lightgreen", "orchid", "gainsboro", "mediumseagreen", "lightgray", "mediumturquoise", "lemonchiffon", "cadetblue", "lightyellow", "lavenderblush", "coral", "purple", "aqua", "whitesmoke", "mediumslateblue", "darkorange", "mediumaquamarine", "darksalmon", "beige", "blueviolet", "azure", "lightsteelblue", "oldlace")
    scales     = ("linear", "log", "symlog", "logit")

    # Fix labels
    texlabels = Dict(
                  "sxx" => L"\sigma_{xx}",
                  "syy" => L"\sigma_{yy}",
                  "szz" => L"\sigma_{zz}",
                  "sxy" => L"\sigma_{xy}",
                  "syz" => L"\sigma_{yz}",
                  "sxz" => L"\sigma_{xz}",
                  "exx" => L"\varepsilon_{xx}",
                  "eyy" => L"\varepsilon_{yy}",
                  "ezz" => L"\varepsilon_{zz}",
                  "exy" => L"\varepsilon_{xy}",
                  "eyz" => L"\varepsilon_{yz}",
                  "exz" => L"\varepsilon_{xz}",
                  "alpha"   => L"\alpha",
                  "beta"    => L"\beta",
                  "gamma"   => L"\gamma",
                  "Gamma"   => L"\Gamma",
                  "delta "  => L"\delta ",
                  "Delta"   => L"\Delta",
                  "epsilon" => L"\epsilon",
                  "zeta"    => L"\zeta",
                  "eta"     => L"\eta",
                  "theta"   => L"\theta",
                  "Theta"   => L"\Theta",
                  "iota"    => L"\iota",
                  "kappa"   => L"\kappa",
                  "lambda"  => L"\lambda",
                  "Lambda"  => L"\Lambda",
                  "mu"      => L"\mu",
                  "nu"      => L"\nu",
                  "omicron" => L"\omicron",
                  "pi"      => L"\pi",
                  "Pi"      => L"\Pi",
                  "rho"     => L"\rho",
                  "sigma"   => L"\sigma",
                  "Sigma"   => L"\Sigma",
                  "tau"     => L"\tau",
                  "upsilon" => L"\upsilon",
                  "Upsilon" => L"\Upsilon",
                  "phi"     => L"\phi",
                  "Phi"     => L"\Phi",
                  "chi"     => L"\chi",
                  "psi"     => L"\psi",
                  "Psi"     => L"\Psi",
                  "omega"   => L"\omega",
                  "Omega"   => L"\Omega",
                 )

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

    # Set axis limits
    xlim!=nothing && plt.xlim(xlim) 
    ylim!=nothing && plt.ylim(ylim) 

    # Print axes labels
    xlabel = get(texlabels, xlabel, xlabel)
    ylabel = get(texlabels, ylabel, ylabel)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Set axis scales
    xscale in scales || error("cplot: xscale should be one of $scales.\n\"$xscale\" was provided.")
    yscale in scales || error("cplot: yscale should be one of $scales.\n\"$yscale\" was provided.")
    plt.xscale(xscale)
    plt.yscale(yscale)

    # Plot curves
    haslegend = false
    for (i, line) in enumerate(data)
        X = get(line, :x, 0)
        Y = get(line, :y, 0)
        lw = get(line, :lw, 0.5)
        ls = get(line, :ls, "-")
        ms = get(line, :ms, 2)
        marker = get(line, :marker, nothing)
        color = get(line, :color, "C$(i-1)")
        label = string(get(line, :label, ""))
        label = get(texlabels, label, label)
        label != "" && (haslegend=true)

        (isa(lw, Number) && lw>0) || error("cplot: lw should be a number greater that zero. '$lw' was provided.")
        (isa(ms, Number) && ms>0) || error("cplot: ms should be a number greater that zero. '$ms' was provided.")
        ls in line_styles || error("cplot: ls should be one of $line_styles.\n\"$ls\" was provided.")
        marker in markers || error("cplot: marker should be one of $markers. \n\"$markers\" was provided")
        color in colors || error("cplot: color should be one of $colors. \n\"$color\" was provided")

        plt.plot(X, Y, marker=marker, ms=ms, color=color, lw=lw, ls=ls, label=label)
    end


    grid && plt.grid( color="lightgrey", ls="dotted", lw=0.3)
    xscale=="linear" && plt.locator_params(axis="x", nbins=xbins)
    yscale=="linear" && plt.locator_params(axis="y", nbins=ybins)

    # plot legend
    if haslegend
        mode = legendexpand ? "expand" : nothing
        ncol = ncol==0 ? length(data) : ncol

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
    if ticksinside
        ax = plt.axes()
        ax.tick_params(axis="x", direction="in")
        ax.tick_params(axis="y", direction="in")
    end
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
