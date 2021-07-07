using LaTeXStrings
export newchart, savechart, showchart, cplot, cplot0


"""
    cplot(data, filename, kwargs...)

Plots a pyplot line chart.

# Arguments

`data` : An array of named tuples with information for each line

`filename` = "" : The chart filename

# Keyword arguments

`xlabel         = "\$x\$"` : x label

`ylabel         = "\$y\$"` : y label

`xbins          = 10` : number of bins in x

`ybins          = 10` : number of bins in y

`grid           = false` : grid

`figsize        = (3.2,2.2)` : size of plot

`legendloc      = "best"` : legend location in plot

`legendexpand   = false` : expanded version

`bbox_to_anchor   = nothing` : legend anchor

`ncol           = 0` : number of columns in horizontal legend

`xlim           = nothing` : x axis limits

`ylim           = nothing` : y axis limits

`xscale         = "linear"` : x axis scale

`yscale         = "linear"` : y axis scale

`xmult          = 1.0` : x data multiplier

`ymult          = 1.0` : y data multiplier

`ticksinside    = true` : put ticks inside plot

`fontsize       = 7` : font size

`legendfontsize = 0` : defaults to fontsize

`labelspacing   = 0.5` : spacing between legend labels

`textlist       = []` : e.g. [ (text="text", pos=(x,y)), ]

`tagpos         = 0.5` : tag position along data

`tagloc         = "top"` :  tag location relative to curve

`tagdist        = 0.005` :  tag distance from curve

`tagfontsize    = 0` : default to fontsize

`tagcolor       = "black"` :  "" uses the data color

`tagalign       = true` : if true the tag is aligned to the curve

`verbosity           = 0` : verbosity level from 0 (silent) to 2 (verbose)

`copypath       = ""` : path to copy ouput file

# Example

```
using Amaru
X1 = collect(1:20)
X2 = collect(1:20)
Y1 = log.(X1)
Y2 = log.(X2)*1.1

cplot([
       (x=X1, y=Y1, marker="o", color="r", label="curve 1")
       (x=X2, y=Y2, marker="s", color="b", label="curve 2")
      ],
      "plot.pdf"
     )
```

Extended example:

```
using Amaru
X1 = collect(1:20)
X2 = collect(1:20)
Y1 = log.(X1)
Y2 = log.(X2)*1.1

cplot([
       (x=X1, y=Y1, lw=0.5, ls="-", ms=2, mfc="w", marker="o", color="r", label="curve 1")
       (x=X2, y=Y2, lw=0.5, ls="-", ms=2, mfc="w", marker="s", color="b", label="curve 2")
      ],
      filename = "plot.pdf",
      xlabel="\$x\$", ylabel="\$y\$", legendloc="best",
      xbins=6, ybins=6, grid=true, figsize=(3,2), legendexpand=false, ncol=0,
      xlim=(0,20), xscale="linear", yscale="log",
      fontsize=7, legendfontsize=0, labelspacing=0.5
     )
```

"""
function cplot(
    data::Array{<:NamedTuple}, 
    fname::String = "";
    filename::String = "",
    xlabel           = L"x",
    ylabel           = L"y",
    xbins            = 10,
    ybins            = 8,
    grid             = true,
    figsize          = (3.2,2.2),
    legendloc        = "best",
    legendexpand     = false,
    bbox_to_anchor   = nothing,
    ncol             = 0,
    xlim             = nothing,
    ylim             = nothing,
    xscale           = "linear",
    yscale           = "linear",
    xmult            = 1.0,
    ymult            = 1.0,
    ticksinside      = true,
    fontsize         = 6.5,
    legendfontsize   = 0,
    labelspacing     = 0.5, 
    annotations      = [],  # e.g. [ (text="C1", pos=(x,y)), ]

    tagpos           = 0.5,
    tagloc           = "top",   # e.g. 1,2,3,4
    tagdist          = 0.005, 
    tagfontsize      = 0,
    tagcolor         = "black",
    tagalign         = true,
    
    y2bins           = 8,
    y2lim            = nothing,
    y2scale          = "linear",
    y2mult           = 1.0,

    copypath         = "",
    verbosity             = false
)

    isempty(filename) && (filename=fname)
    verbosity = clamp(verbosity, 0,2)

    verbosity>0 && headline("Chart plotting")
    verbosity>0 && message("generating plot $(strip(xlabel,'$')) vs $(strip(ylabel,'$'))")

    if verbosity>0
        hint("Optional arguments:", level=2)
        options = "xlabel, ylabel, legendloc, xbins, ybins, grid, figsize, legendexpand, ncol,
                xlim, ylim, xscale, yscale, fontsize, ticksinside, legendfontsize, labelspacing"
        hint(options, level=3)

        hint("Arguments and optional arguments per curve:", level=2)
        options = "x, y, color, ls, lw, marker, ms, mfc, label, tag, tagpos, tagloc, tagdist, tagfontsize, tagcolor, tagalign, at_y2"
        hint(options, level=3)
    end

    @eval import PyPlot:plt, matplotlib, figure, ioff, gcf

    line_styles = ("-", "--", "-.", ":", "", " ", "None", nothing)
    markers     = (".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "*", "h", "H", "+", "x", "D", "d", "|", "_", "P", "X", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, "None", nothing, " ", "")
    colors      = ("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
                   "c", "b", "w", "g", "y", "k", "r", "m",
                   "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
                   "indigo", "gold", "hotpink", "firebrick", "indianred", "yellow", "mistyrose", "darkolivegreen", "olive", "darkseagreen", "pink", "tomato", "lightcoral", "orangered", "navajowhite", "lime", "palegreen", "darkslategrey", "greenyellow", "burlywood", "seashell", "mediumspringgreen", "fuchsia", "papayawhip", "blanchedalmond", "chartreuse", "dimgray", "black", "peachpuff", "springgreen", "aquamarine", "white", "orange", "lightsalmon", "darkslategray", "brown", "ivory", "dodgerblue", "per", "lawngreen", "chocolate", "crimson", "forestgreen", "darkgrey", "lightseagreen", "cyan", "mintcream", "silver", "antiquewhite", "mediumorchid", "skyblue", "gray", "darkturquoise", "goldenrod", "darkgreen", "floralwhite", "darkviolet", "darkgray", "moccasin", "saddlebrown", "grey", "darkslateblue", "lightskyblue", "lightpink", "mediumvioletred", "slategrey", "red", "deeppink", "limegreen", "darkmagenta", "palegoldenrod", "plum", "turquoise", "lightgrey", "lightgoldenrodyellow", "darkgoldenrod", "lavender", "maroon", "yellowgreen", "sandybrown", "thistle", "violet", "navy", "magenta", "dimgrey", "tan", "rosybrown", "olivedrab", "blue", "lightblue", "ghostwhite", "honeydew", "cornflowerblue", "slateblue", "linen", "darkblue", "powderblue", "seagreen", "darkkhaki", "snow", "sienna", "mediumblue", "royalblue", "lightcyan", "green", "mediumpurple", "midnightblue", "cornsilk", "paleturquoise", "bisque", "slategray", "darkcyan", "khaki", "wheat", "teal", "darkorchid", "salmon", "deepskyblue", "rebeccapurple", "darkred", "steelblue", "palevioletred", "lightslategray", "aliceblue", "lightslategrey", "lightgreen", "orchid", "gainsboro", "mediumseagreen", "lightgray", "mediumturquoise", "lemonchiffon", "cadetblue", "lightyellow", "lavenderblush", "coral", "purple", "aqua", "whitesmoke", "mediumslateblue", "darkorange", "mediumaquamarine", "darksalmon", "beige", "blueviolet", "azure", "lightsteelblue", "oldlace")
    scales      = ("linear", "log", "symlog", "logit")

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

    #=
    import matplotlib.font_manager as font_manager

    font_dirs  = ['/my/custom/font/dir', ]
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    font_list  = font_manager.createFontList(font_files)
    font_manager.fontManager.ttflist.extend(font_list)

    mpl.rcParams['font.family'] = 'My Custom Font'
    =#

    # Configure plot
    legendfontsize==0 && (legendfontsize=fontsize)
    tagfontsize==0 && (tagfontsize=fontsize)

    if filename!=""
        plt.rc("font", family="STIXGeneral", size=fontsize)
        plt.rc("mathtext", fontset="cm")
        plt.rc("lines", scale_dashes=true)

        plt.close("all")
        plt.ioff()
        plt.rc("xtick", labelsize=fontsize)
        plt.rc("ytick", labelsize=fontsize)
        plt.rc("lines", lw=0.7)
        plt.rc("lines", markersize=1.5)
        plt.rc("axes" , linewidth=0.5)
        plt.rc("figure", figsize=figsize)
        plt.rc("legend", fontsize=legendfontsize)
    else
        plt.rc("font", family="STIXGeneral", size=fontsize+3)
        plt.rc("mathtext", fontset="cm")
        plt.ion()
    end

    # Newx axes
    # ax = plt.axes()

    # Set axis limits
    xlim!==nothing && plt.xlim(xlim) 
    ylim!==nothing && plt.ylim(ylim) 

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
        X = get(line, :x, 0).*xmult
        Y = get(line, :y, 0).*ymult
        length(X) == length(Y) || error("cplot: (line $i) x and y lengths do not match")
        lw     = get(line, :lw, 0.5)
        ls     = get(line, :ls, "-")
        color  = get(line, :color, "C$(i-1)")
        marker = get(line, :marker, nothing)
        ms     = get(line, :ms, 2)
        mfc    = get(line, :mfc, color)
        label  = string(get(line, :label, ""))
        label  = get(texlabels, label, label)
        label != "" && (haslegend=true)

        (isa(lw, Number) && lw>0) || error("cplot: (line $i) lw should be a number greater that zero. '$lw' was provided.")
        (isa(ms, Number) && ms>0) || error("cplot: (line $i) ms should be a number greater that zero. '$ms' was provided.")
        ls in line_styles || error("cplot: (line $i) ls should be one of $line_styles.\n\"$ls\" was provided.")
        marker in markers || error("cplot: (line $i) marker should be one of $markers. \n\"$markers\" was provided")
        #(color isa Tuple || color in colors) || error("cplot: color should be one of $colors. \n\"$color\" was provided")

        plt.plot(X, Y, marker=marker, ms=ms, mfc=mfc, mew=0.5, color=color, lw=lw, ls=ls, label=label)
    end


    grid && plt.grid(color="lightgrey", which="both", ls="dotted", lw=0.3)
    xscale=="linear" && plt.locator_params(axis="x", nbins=xbins)
    yscale=="linear" && plt.locator_params(axis="y", nbins=ybins)

    # plot legend
    if haslegend
        mode = legendexpand ? "expand" : nothing
        if ncol==0
            if legendloc in ("top", "right", "bottom")
                ncol = length(data)
            else
                ncol = 1
            end
        end

        if legendloc=="top"
            leg = plt.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.01, 1.04, 0.2), edgecolor="k", ncol=ncol, mode=mode)
        elseif legendloc=="right"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1), edgecolor="k")
        elseif legendloc=="bottom"
            leg = plt.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.02, 1.04, -0.2), edgecolor="k", ncol=ncol, mode=mode)
        else
            leg = plt.legend(loc=legendloc, bbox_to_anchor=bbox_to_anchor, edgecolor="k", labelspacing=labelspacing, ncol=ncol)
        end

        frame = leg.get_frame()
        frame.set_linewidth(0.5)
    end

    # Tick parameters
    ax = plt.gca()
    ax.xaxis.set_tick_params(width=0.3)
    ax.yaxis.set_tick_params(width=0.3)
    ax.xaxis.set_tick_params(size=2.5)
    ax.yaxis.set_tick_params(size=2.5)
    if ticksinside
        ax.tick_params(which="minor", axis="x", direction="in")
        ax.tick_params(which="minor", axis="y", direction="in")
        ax.tick_params(which="major", axis="x", direction="in")
        ax.tick_params(which="major", axis="y", direction="in")
    end

    # print text
    if length(annotations)>0
        for line in annotations
            x, y     = get(line, :pos, (0.5, 0.5))
            fontsize = get(line, :fontsize, tagfontsize)
            text     = string(get(line, :text, ""))
            arrowcoord = get(line, :arrowcoord, nothing)
            
            if arrowcoord===nothing
                plt.text(x, y, text, fontsize=fontsize, transform=ax.transAxes, ha="center", va="center")
            else
                plt.annotate(text, 
                    xytext     = (x,y), 
                    textcoords = "figure fraction",
                    xy         = (arrowcoord[1]*xmult, arrowcoord[2]*ymult),
                    arrowprops = Dict(
                        "lw"              => 0.3,
                        "arrowstyle"      => "->",
                        "facecolor"       => "black",
                        "connectionstyle" => "angle,angleA=0,angleB=90,rad=3"
                    ),

                )
            end
            
        end
    end

    # line tags
    # from Axes to Data coords:  (ax.transAxes+ax.transData.inverted()).transform([x,y])
    # from Data to Axes coords:  (ax.transData+ax.transAxes.inverted()).transform([x,y])
    for (i, line) in enumerate(data)
        tag = get(line, :tag, "")
        tag != "" || continue

        X     = get(line, :x, 0).*xmult
        Y     = get(line, :y, 0).*ymult
        pos   = get(line, :tagpos, tagpos)
        loc   = get(line, :tagloc, tagloc)
        dist  = get(line, :tagdist, tagdist)
        align = get(line, :tagalign, tagalign)
        color = get(line, :tagcolor, tagcolor)
        
        if color==""
            color = get(line, :color, "C$(i-1)")
            color = 0.85.*matplotlib.colors.to_rgb(color)
        end

        XY    = (ax.transData+ax.transAxes.inverted()).transform([X Y])
        X     = XY[:,1]
        Y     = XY[:,2]

        x = y = 0.0
        dx = dy = 0.0
        α = 0.0
        # @show X
        if pos isa Array
            n = length(X)
            m = length(pos)
            found = false # j
            for j=2:m
                A1x, A1y = pos[j-1]
                A2x, A2y = pos[j] .+ 1e-3

                for i=2:n
                    B1x, B1y = X[i-1], Y[i-1]
                    B2x, B2y = X[i], Y[i]

                    max(A1x, A2x) < min(B1x, B2x)  && continue
                    min(A1x, A2x) > max(B1x, B2x)  && continue

                    M = [ A2x-A1x B1x-B2x
                          A2y-A1y B1y-B2y ]
                    r, s = inv(M)*[ B1x-A1x, B1y-A1y ]
                    x, y = (A1x + r*(A2x-A1x), A1y + r*(A2y-A1y))

                    if A1x <= x <= A2x && B1x <= x <= B2x
                        α = atand(Y[i]-Y[i-1], X[i]-X[i-1])
                        dx = dist*abs(cosd(α))
                        dy = dist*abs(sind(α))
                        found = true
                        break
                    end
                end

                found && break
            end
        else # pos is a scalar
            len = 0.0
            for i=2:length(X)
                len += √((X[i]-X[i-1])^2 + (Y[i]-Y[i-1])^2)
            end
            lpos = pos*len
            len = 0.0

            for i=2:length(X)
                dlen = √((X[i]-X[i-1])^2 + (Y[i]-Y[i-1])^2)
                len += dlen
                if len>=lpos
                    x = X[i] - (len-lpos)/dlen*(X[i]-X[i-1])
                    y = Y[i] - (len-lpos)/dlen*(Y[i]-Y[i-1])
                    α = atand(Y[i]-Y[i-1], X[i]-X[i-1])
                    dx = dist*abs(cosd(α))
                    dy = dist*abs(sind(α))
                    break
                end
            end
        end


        # Default location "top"
        if loc in ("top", "t", 't')
            if 0 < α <= 90 || -180 < α <= -90
                ha = "right"; va = "bottom"
                dx, dy = -dy, dx
            elseif 90 < α <= 180 || -90 < α <= 0
                ha = "left"; va = "bottom"
                dx, dy = dy, dx
            end
        else
            if 0 < α <= 90 || -180 < α <= -90
                ha = "left"; va = "top"
                dx, dy = dy, -dx*5
            elseif 90 < α <= 180 || -90 < α <= 0
                ha = "right"; va = "top"
                dx, dy = -dy, -dx*5
            end
        end

        xpos, ypos = x+dx, y+dy

        if align
            align && (ha="center")
            
            # println()
            # @s α
            90<α<=180 && (α=α-180)
            -180<α<=-90 && (α=α+180)
            # @s α
            angle = ax.transAxes.transform_angles([α], [1.0 1.0])[1]
        else
            angle = 0.0
        end


        # plt.text(xpos, ypos, tag, ha=ha, va=va, color=color, fontsize=tagfontsize, transform=ax.transAxes)
        plt.text(
            xpos, ypos, tag, ha=ha, va=va, color=color, 
            # backgroundcolor="white",
            fontsize=tagfontsize, 
            rotation=angle, rotation_mode="anchor",
            transform=ax.transAxes
        )
    end

    # second axis
    # ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    # ax2.plot(t, data2, color=color)
    # ax2.set_ylim([0, 2.8])
    # ax2.set_ylabel("Average stress [MPa]")  # we already handled the x-label with ax1
    # ax2.tick_params(axis='y', labelcolor=color)
     
    # show or save plot
    if filename=="" 
        plt.fignum_exists(ax.figure.number) || plt.show()
    else
        _, format = splitext(filename)
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.01, format=format[2:end])
        verbosity>0 && info("file $filename saved")
        plt.close("all")

        if copypath!=""
            if isdir(copypath)
                copyfile = joinpath(copypath, basename(filename))
            else
                copyfile = copypath
            end
            try
                cp(filename, copyfile, force=true)
                verbosity>0 && info("file $copyfile saved")
            catch err
                notify("cplot: $filename could not be copied to $copypath")
            end
        end
    end

    return @eval gcf()

end