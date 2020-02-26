
function newfig(; xlabel="\$x\$", ylabel="\$y\$", lw=0.7, ms=2, legendloc="best",
               xbins=6, ybins=6, grid=false, figsize=(3,2), legendexpand=false, ncol=0,
               xmin=NaN, xmax=NaN, ymin=NaN, ymax=NaN,
               fontsize=7, legendfontsize=0, labelspacing=0.5)

    @eval import PyPlot:plt

    # Configure plot
    plt.close("all")

    plt.rc("font", family="STIXGeneral", size=fontsize)
    plt.rc("mathtext", fontset="cm")
    plt.rc("lines", scale_dashes=true)

    plt.ioff()
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

    return nothing
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
