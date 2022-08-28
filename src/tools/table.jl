# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DataTable, datatable, push!, save, loadtable, randtable
export compress, resize, filter, cut!, clamp!, denoise!
export getheader


# DataTable object
const KeyType = Union{Symbol,AbstractString}

mutable struct DataTable
    columns::Vector{AbstractVector}
    colidx ::OrderedDict{String,Int} # Data index
    header ::Array{String,1}
    name   ::String
    function DataTable(header::Array, columns::Array{<:AbstractVector,1}=Array{Any,1}[]; name="")
        @assert length(header) == length(unique(header))
        if length(columns)==0
            columns = AbstractVector[ [] for i in 1:length(header) ]
        end

        colidx = OrderedDict{String,Int}( string(key)=>i for (i,key) in enumerate(header) )

        return new(columns, colidx, string.(header), name)
    end
end

getcolumns(table::DataTable) = getfield(table, :columns)
getheader(table::DataTable) = getfield(table, :header)
getcolidx(table::DataTable) = getfield(table, :colidx)
getname(table::DataTable) = getfield(table, :name)


function DataTable(; pairs...)

    columns = AbstractVector[]
    header  = KeyType[]
    for (key, val) in pairs
        push!(header, string(key))
        push!(columns, copy(val))
    end

    return DataTable(header, columns)
end


function DataTable(header::Array, matrix::Array{T,2} where T; name="")
    nkeys = length(header)
    ncols = size(matrix,2)
    nkeys != ncols && error("DataTable: header and number of data columns do not match")
    types = [ typeof(matrix[1,i]) for i in 1:ncols ]
    columns = AbstractVector[ convert(Array{types[i],1}, matrix[:,i]) for i in 1:ncols ]
    return DataTable(header, columns, name=name)
end


function ncols(table::DataTable)
    return length(getcolumns(table))
end


function nrows(table::DataTable)
    columns = getcolumns(table)
    ncols   = length(columns)
    ncols == 0 && return 0
    return length(columns[1])
end


function Base.size(table::DataTable)
    columns = getcolumns(table)
    ncols   = length(columns)
    ncols == 0 && return (0,0)
    return (length(columns[1]), ncols)
end


function Base.getproperty(table::DataTable, sym::Symbol)
    return getindex(table, sym)
end


function Base.setproperty!(table::DataTable, sym::Symbol, val)
    key = string(sym)
    setindex!(table, val, key)
end


function Base.push!(table::DataTable, row::Array{T,1} where T)
    nr, nc = size(table)
    @assert nc==length(row)
    columns = getcolumns(table)

    for (i,v) in enumerate(row)
        push!(columns[i], v)
    end
end


function Base.keys(table::DataTable)
    return keys(getcolidx(table))
end


function Base.push!(table::DataTable, dict::AbstractDict)
    columns = getcolumns(table)
    colidx  = getcolidx(table)
    header  = getheader(table)

    if length(columns)==0
        for (k,v) in dict
            key = string(k)
            push!(header, key)
            push!(columns, [v])
            colidx[key] = length(header)
        end
    else
        nrows = length(columns[1])
        for (k,v) in dict
            key = string(k)
            # Add data
            idx = get(colidx, key, 0)
            if idx==0
                # add new column
                new_col = zeros(nrows)
                push!(new_col, v)
                push!(columns, new_col)
                colidx[string(k)] = length(columns)
                push!(header, key)
            else
                push!(columns[idx], v)
            end
        end

        # Add zero for missing values
        for col in columns
            if length(col)==nrows
                push!(col, 0.0)
            end
        end
    end
end


function Base.setindex!(table::DataTable, column::AbstractVector, key::KeyType)
    columns  = getcolumns(table)
    colidx = getcolidx(table)
    header   = getheader(table)

    if length(columns)>0
        n1 = length(columns[1])
        n2 = length(column)
        n1>0 && n1!=n2 && error("setindex! : vector length ($n2) for data ($key) is incompatible with DataTable rows ($n1)")
    end

    key = string(key)
    if haskey(colidx, key)
        idx = colidx[key]
        columns[idx] = column
    else
        push!(columns, column)
        push!(header, key)
        colidx[key] = length(columns)
    end
    return column
end


function Base.getindex(table::DataTable, key::KeyType)
    key    = string(key)
    colidx = getcolidx(table)
    haskey(colidx, key) || error("getindex: key $(repr(key)) not found in DataTable. Available keys are $(keys(table))")

    idx  = colidx[key]
    return getcolumns(table)[idx]
end

function Base.getindex(table::DataTable, keys::Array{<:KeyType,1})
    columns = [ table[string(key)] for key in keys ]
    return DataTable(keys, columns)
end

function Base.getindex(table::DataTable, rowindex::Int)
    columns = getcolumns(table)
    cols = [ [ columns[i][rowindex] ] for i in 1:length(columns) ]

    return DataTable(getheader(table), cols)
end

function Base.getindex(table::DataTable, idxs::Union{Colon,OrdinalRange{Int,Int},Array{Int,1},BitArray{1}})
    columns = getcolumns(table)
    cols = [ columns[i][idxs] for i in 1:length(columns) ]

    return DataTable(getheader(table), cols)
end


function Base.getindex(table::DataTable, rows, col::KeyType)
    return table[rows][col][1]
end

function Base.getindex(table::DataTable, rows, cols::Array{<:KeyType,1})
    return table[rows][cols]
end




function Base.Array(table::DataTable)
    return hcat(getcolumns(table)...)    
end


# function Base.getindex(table::DataTable, rowindex::Int, colidx::Int)
#     return getcolumns(table)[colidx][rowindex]
# end

function Base.getindex(table::DataTable, rowindex::Union{Int,Colon,OrdinalRange{Int,Int},Array{Int,1}}, colidx::Int)
    return getcolumns(table)[colidx][rowindex]
end


function Base.getindex(table::DataTable, rowindex::Int, colon::Colon)
    return [ column[rowindex] for column in getcolumns(table)]
end

# function Base.getindex(table::DataTable, rows, cols::Union{Colon,OrdinalRange{Int,Int},Array{Int,1}})
#     return return getcolumns(table)[colidx][rows]
# end


function Base.lastindex(table::DataTable)
    nr, nc = size(table)
    nc==0 && error("lastindex: use of 'end' in an empty DataTable")
    return nr
end


function Base.lastindex(table::DataTable, idx::Int)
    nr, nc = size(table)
    nc==0 && error("lastindex: use of 'end' in an empty DataTable")
    if idx==1
        return nr
    else
        return nc
    end
end


sprintf(fmt, args...) = @eval @sprintf($fmt, $(args...))



function Base.filter(table::DataTable, expr::Expr)
    fields = get_vars(expr)
    vars   = Dict{Symbol, Float64}()
    nr     = nrows(table)
    idxs   = falses(nr)
    for i in 1:nr
        for field in fields
            vars[field] = table[field][i]
        end
        idxs[i] = eval_arith_expr(expr; vars...)
    end
    return table[idxs]
end


function Base.sort!(table::DataTable, options::NamedTuple...)
    n, m = size(table)
    idx = collect(1:n)
    
    for opt in options
        field = get(opt, :field, "")
        rev   = get(opt, :rev, false)

        col  = table[field][idx]
        idx_ = sortperm(col, rev=rev)
        idx  = idx[idx_]
    end

    cols = getcolumns(table)
    for i in 1:m
        cols[i] = cols[i][idx]
    end

    return table
end


function compress!(table::DataTable, n::Int)
    columns = getcolumns(table)
    nr, nc = size(table)

    nr<=n && return table[:]

    factor = (nr-1)/(n-1)
    idxs   = [ round(Int, 1+factor*(i-1)) for i in 1:n ]

    for i in 1:length(columns)
        columns[i] = columns[i][idxs]
    end
    
    return table
end


function resize(table::DataTable, n::Int=0; ratio=1.0)
    header = getheader(table)
    nr     = nrows(table)
    
    if n==0
        ratio > 0.0 || error("resize: ratio should be greater than zero")
        n = max(2, round(Int, nr*ratio))
    end
    n > 0 || error("resize: number of rows should be greater than zero")
    nr >= 4 || error("resize: Table object should contain at least 4 rows")

    ns = nr - 1                # number of spacings
    nb = floor(Int64, ns/3)    # number of bezier curves
    ds = 1.0 / nb              # spacing between curves

    cols = Array{Any,1}[]
    for (j,field) in enumerate(header)
        U = table[field]
        V = zeros(n)

        for (i,s) in enumerate(range(0.0, 1.0, length=n))
            # find index of Bezier and local coordinate t
            ib = floor(Int64, s/ds) + 1   # index of Bezier
            ib > nb && (ib = nb)          # fix index if s ~= 1+eps
            s0 = (ib-1) * ds              # s @ left point
            t  = (s - s0) / ds            # local t for current Bezier
            t > 1.0 && (t = 1.0)          # clean rubbish. e.g. 1.00000000002
        
            # collect control points
            k = 1 + (ib-1) * 3            # position of first point of bezier
            
            P1 = U[k  ]
            P2 = U[k+1]
            P3 = U[k+2]
            P4 = U[k+3]

            # control points
            Q1 =         P1
            Q2 = (-5.0 * P1 + 18.0 * P2 -  9.0 * P3 + 2.0 * P4) / 6.0
            Q3 = ( 2.0 * P1 -  9.0 * P2 + 18.0 * P3 - 5.0 * P4) / 6.0
            Q4 =                                            P4

            a =       Q4 - 3.0 * Q3 + 3.0 * Q2 - Q1
            b = 3.0 * Q3 - 6.0 * Q2 + 3.0 * Q1
            c = 3.0 * Q2 - 3.0 * Q1
            d =       Q1

            V[i] = a*t*t*t + b*t*t + c*t + d
        end

        push!(cols, V)
    end
    
    return DataTable(header, cols)
end


function cut!(table::DataTable, field, value=0.0; after=false)
    V   = table[field] .- value
    idx = 0
    for i in 2:length(V)
        if V[i-1]*V[i] <= 0
            idx = i
            break
        end
    end

    if idx>0
        α = -V[idx-1]/(V[idx]-V[idx-1]) 
        header = getheader(table)
        for field in header
            W      = table[field]
            W[idx] = W[idx-1] + α*(W[idx]-W[idx-1])
        end
        rng = 1:idx
        after && (rng=1:length(V))

        columns = getcolumns(table)
        for (i,col) in enumerate(columns)
            columns[i] = col[rng]
        end
    end

    return table
end


function clamp!(table::DataTable, field, lo, hi)
    clamp!(table[field], lo, hi)
    return table
end

export smooth!
function smooth!(table::DataTable, fieldx, fieldy=nothing; knots=[0.0, 1.0])
    nr     = nrows(table)

    if fieldy === nothing
        fieldy = fieldx
        X = collect(range(0, 1, length=nr))
        Y = table[fieldx]
    else 
        X = table[fieldx]
        Y = table[fieldy]
    end

    Idxs = split(X, X[1] .+ knots.*(X[end]-X[1]))
    # @show X
    # @show Idxs

    for i in 1:length(Idxs)
        Xi = X[Idxs[i]]
        Yi = Y[Idxs[i]]
        M = Float64[ ones(length(Xi)) Xi Xi.^3 ]
        a, b, d = inv(M'*M)*M'*Yi
        Y[Idxs[i]] .= a .+ b*Xi .+ d*Xi.^3
    end

    table[fieldy] = Y
end

function denoise!(table::DataTable, fieldx, fieldy=nothing; noise=0.05, npatch=4)
    header = getheader(table)
    nr     = nrows(table)

    if fieldy === nothing
        X = range(0,1,length=nr)
        Y = table[fieldx]
    else
        X = table[fieldx]
        Y = table[fieldy]
    end

    # Regression along patches
    M  = ones(npatch,2)
    A  = zeros(2)
    ΔY = [ Float64[] for i in 1:nr ]

    for i in 1:nr-npatch+1        
        rng = i:i+npatch-1
        Xp = X[rng]
        Yp = Y[rng]

        # Linear regression
        M[:,2] .= Xp
        A   = pinv(M)*Yp
        ΔYp = abs.(Yp .- M*A)

        for (j,k) in enumerate(rng)
            push!(ΔY[k], ΔYp[j])
        end
    end

    # Get indexes to keep
    idxs = minimum.(ΔY) .<= noise*(maximum(Y)-minimum(Y))
    # idxs[1:npatch] .= 1
    # idxs[end-npatch+1:end] .= 1

    # newtable = table[:]
    # for i in (1:nr)[idxs]
    #     j = findprev(!iszero, idxs, i-1)
    #     k = findnext(!iszero, idxs, i+1)
    #     r = (X[i]-X[j])/(X[k]-X[j])

    #     for fieldx in header
    #         V = newtable[fieldx]
    #         V[i] = V[j] + r*(V[k]-V[j])
    #     end
    # end

    # Linear interpolation of dropped points
    cols = Array{Any,1}[]
    for column in getcolumns(table)
        # V = copy(column)
        V = column
        for i in (1:nr)[.!idxs]
            j = findprev(!iszero, idxs, i-1)
            k = findnext(!iszero, idxs, i+1)
            if j===nothing
                j = k
                k = findnext(!iszero, idxs, j+1)
            end
            if k===nothing
                k = j
                j = findprev(!iszero, idxs, k-1)
            end
            r = (X[i]-X[j])/(X[k]-X[j])
            V[i] = V[j] + r*(V[k]-V[j])
        end
        push!(cols, V)
    end

    return table
    # return DataTable(header, cols)
end


# TODO: Improve column width for string items
function save(table::DataTable, filename::String; report=false, digits::Array=[])
    suitable_formats = (".dat", ".table", ".tex")

    basename, format = splitext(filename)
    format in suitable_formats || error("DataTable: cannot save in \"$format\" format. Suitable formats $suitable_formats.")

    local f::IOStream
    try
        f  = open(filename, "w")
    catch err
        warn("DataTable: File $filename could not be opened for writing.")
        return
    end

    columns = getcolumns(table)
    colidx  = getcolidx(table)
    header  = getheader(table)
    nr, nc  = size(table)

    if format in (".dat", ".table")
        for (i,key) in enumerate(header)
            @printf(f, "%12s", key)
            print(f, i!=nc ? "\t" : "\n")
        end

        # print values
        for i in 1:nr
            for j in 1:nc
                item = columns[j][i]
                if typeof(item)<:AbstractFloat
                    @printf(f, "%12.5e", item)
                elseif typeof(item)<:Integer
                    @printf(f, "%12d", item)
                else
                    @printf(f, "%12s", item)
                end
                print(f, j!=nc ? "\t" : "\n")
            end
        end

        report && printstyled("  file $filename written\n", color=:cyan)
    end

    if format==".tex"
        # widths calculation
        widths = length.(header)
        types  = eltype.(columns)

        if length(digits)==0
            digits = repeat([4], nc)
        end
        @assert length(digits)==nc

        for (i,col) in enumerate(columns)
            etype = types[i]
            if etype<:AbstractFloat
                widths[i] = max(widths[i], 12)
            elseif etype<:Integer
                widths[i] = max(widths[i], 6)
            elseif etype<:AbstractString
                widths[i] = max(widths[i], maximum(length.(col)))
            else
                widths[i] = max(widths[i], maximum(length.(string.(col))))
            end
        end

        # printing header
        level = 1
        indent = "    "
        println(f, indent^level, raw"\begin{tabular}{", "c"^nc, "}" )
        level = 2
        println(f, indent^level, raw"\toprule")
        print(f, indent^level)
        for (i,key) in enumerate(header)
            etype = types[i]
            width = widths[i]
            if etype<:Real
                print(f, lpad(key, width))
            else
                print(f, rpad(key, width))
            end
            i<nc && print(f, " & ")
        end
        println(f, raw" \\\\")

        # printing body
        println(f, indent^level, raw"\hline")
        for i in 1:nr
            print(f, indent^level)
            for j in 1:nc
                etype = types[j]
                item = columns[j][i]
                width = widths[j]
                if etype<:AbstractFloat
                    #item = @sprintf("%12.3f", item)
                    dig = digits[j]
                    if isnan(item)
                        item = "-"
                    else
                        item = sprintf("%$width.$(dig)f", item)
                    end
                    print(f, lpad(string(item), width))
                elseif etype<:Integer
                    item = @sprintf("%6d", item)
                    print(f, lpad(item, width))
                elseif etype<:AbstractString
                    print(f,rpad(item, width))
                else
                    str = string(item)
                    print(f, rpad(item, width))
                end
                j<nc && print(f, " & ")
            end
            println(f, raw" \\\\")
        end
        println(f, indent^level, raw"\bottomrule")

        # printing ending
        level = 1
        println(f, indent^level, raw"\end{tabular}")
    end

    close(f)
    return nothing
end


function DataTable(filename::String, delim='\t')
    basename, format = splitext(filename)
    formats = (".dat", ".table")
    format in formats || error("DataTable: cannot read \"$format\". Suitable formats are $formats")

    if format in (".dat", ".table")
        matrix, headstr = readdlm(filename, delim, header=true, use_mmap=false)
        header = vec(strip.(headstr))
        table = DataTable(header, matrix)
        return table
    end
end


# Functions for backwards compatibility
loadtable(filename::String, delim='\t') = DataTable(filename, delim)


# TODO: Improve display. Include column datatype
function Base.show(io::IO, table::DataTable)
    columns = getcolumns(table)
    colidx = getcolidx(table)

    println(io)
    if length(columns)==0
        print(io, "DataTable()")
        return
    end

    nr, nc = size(table)

    if nr==0
        print(io, "DataTable()")
        return
    end

    header = keys(colidx)
    types  = typeof.(getindex.(columns,1))

    hwidths = length.(header)
    widths  = zeros(Int, length(header))
    useformat   = falses(length(header))
    shortformat = falses(length(header))

    for (i,col) in enumerate(columns)
        etype = types[i]
        widths[i] = maximum(length.(string.(col)))
        if etype<:AbstractFloat
            if widths[i] >= 11
                widths[i] = 11
                useformat[i] = true
            end
        end
        widths[i] = max(widths[i], hwidths[i])
    end

    total = sum(widths) + (length(header)+1)*3
    if total>displaysize(stdout)[2]
        for (i,col) in enumerate(columns)
            etype = types[i]
            if etype<:AbstractFloat && useformat[i] && hwidths[i]<=8
                shortformat[i] = true
                widths[i] = max(8, hwidths[i])
            end
        end
    end

    print(io, " │ ")
    for (i,key) in enumerate(header)
        etype = types[i]
        width = widths[i]
        if etype<:Real
            print(io, lpad(key, width))
        else
            print(io, rpad(key, width))
        end
        print(io, " │ ")
    end
    println(io)

    visible_rows = 30
    half_vrows = div(visible_rows,2)

    # print values
    for i in 1:nr
        if i>half_vrows && nr-i>=half_vrows
            i==half_vrows+1 && println(io, " ⋮")
            continue
        end

        print(io, " │ ")
        for j in 1:nc
            etype = types[j]
            item = columns[j][i]
            if etype<:AbstractFloat
                if useformat[j]
                    if shortformat[j]
                        item = @sprintf("%8.1e", item)
                    else
                        item = @sprintf("%11.4e", item)
                    end
                end
                print(io, lpad(item, widths[j]))
            else
                print(io, rpad(item, widths[j]))
            end
            print(io, " │ ")
        end
        i<nr && println(io)
    end

end


randtable() = DataTable(["x","y","z"], [0:10 rand().*(sin.(0:10).+(0:10)) rand().*(cos.(0:10).+(0:10)) ])
