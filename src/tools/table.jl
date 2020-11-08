# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DataTable, DataBook, push!, save, loadtable, loadbook, randtable, compress, resize, filter


# DataTable object
const KeyType = Union{Symbol,AbstractString}
const ColType = Array{T,1} where T

mutable struct DataTable
    columns  :: Array{ColType,1}
    colindex :: OrderedDict{String,Int} # Data index
    header   :: Array{String,1}
    name     :: String
    function DataTable(header::Array)
        this = new()
        header = vec(header)
        this.columns  = [ [] for s in header ]
        this.colindex = OrderedDict( string(key)=>i for (i,key) in enumerate(header) )
        this.header   = string.(header)
        return this
    end
end


function DataTable(header::Array, columns::Array{<:ColType,1})
    this      = DataTable(header)
    nfields   = length(header)
    ncols     = length(columns)
    nfields  != ncols && error("DataTable: header and number of data columns do not match")
    this.columns = deepcopy(columns)
    return this
end


function DataTable(; kwargs...)

    columns = ColType[]
    header  = KeyType[]
    for (key, val) in kwargs
        push!(header, string(key))
        push!(columns, copy(val))
    end

    return DataTable(header, columns)
end

function DataTable(header::Array, matrix::Array{T,2} where T)
    this   = DataTable(header)
    nkeys  = length(header)
    ncols  = size(matrix,2)
    nkeys != ncols && error("DataTable: header and number of data columns do not match")
    types = [ typeof(matrix[1,i]) for i=1:ncols ]
    this.columns = [ convert(Array{types[i],1}, matrix[:,i]) for i=1:ncols ]
    return this
end


import Base.push!
function push!(table::DataTable, row::Array{T,1} where T)
    @assert length(table.colindex)==length(row)

    if length(table.columns[1])==0
        table.columns = [ typeof(v)[v] for v in row  ]
    else
        for (i,val) in enumerate(row)
            push!(table.columns[i], val)
        end
    end
end


function Base.keys(table::DataTable)
    return keys(table.colindex)
end

function Base.push!(table::DataTable, dict::AbstractDict)
    if length(table.columns)==0
        table.columns  = [ typeof(v)[v] for (k,v) in dict ]
        table.colindex = OrderedDict( string(key)=>i for (i,key) in enumerate(keys(dict)) )
        table.header   = string.(keys(dict))
    else
        nrows = length(table.columns[1])
        for (k,v) in dict
            # Add data
            colindex = get(table.colindex, string(k), 0)
            if colindex==0
                # add new column
                new_col = zeros(nrows)
                push!(new_col, v)
                push!(table.columns, new_col)
                table.colindex[string(k)] = length(table.columns)
                push!(table.header, string(k))
            else
                push!(table.columns[colindex], v)
            end
        end

        # Add zero for missing values if any
        for col in table.columns
            if length(col)==nrows
                push!(col, 0.0)
            end
        end
    end
end


function Base.setindex!(table::DataTable, column::ColType, key::KeyType)
    col = column

    if length(table.columns)>0
        n1 = length(table.columns[1])
        n2 = length(column)
        n1!=n2 && error("setindex! : length ($n2) for data ($key) is incompatible with DataTable rows ($n1)")
    end

    key = string(key)
    if haskey(table.colindex, key)
        idx = table.colindex[key]
        table.columns[idx] = column
    else
        push!(table.columns, column)
        push!(table.header, key)
        table.colindex[key] = length(table.columns)
    end
    return column
end


function Base.getindex(table::DataTable, key::KeyType)
    string(key) in table.header || error("getindex: key ($(repr(key))) not found in DataTable header. Available keys: $(keys(table))")
    return table.columns[table.colindex[string(key)]]
end

function Base.getindex(table::DataTable, keys::Array{<:KeyType,1})
    columns = [ table[string(key)] for key in keys ]
    subtable = DataTable(keys, columns)
    return subtable
end

function Base.getindex(table::DataTable, rowindex::Int)
    subtable = DataTable(table.header)
    for i in 1:length(table.columns)
        push!(subtable.columns[i], table.columns[i][rowindex])
    end
    return subtable
end

function Base.getindex(table::DataTable, idxs::Union{Colon,OrdinalRange{Int,Int},Array{Int,1}})
    subtable = DataTable(table.header)
    for i in 1:length(table.columns)
        subtable.columns[i] = table.columns[i][idxs]
    end
    return subtable
end


function Base.getindex(table::DataTable, rows::Union{Colon,OrdinalRange{Int,Int},Int,Array{Int,1}}, keys::KeyType)
    if typeof(keys) <: KeyType
        keys = [ keys ]
    end
    return table[keys][rows]
end

function Base.getindex(table::DataTable, rows, cols)
    return table[cols][rows]
end

function Base.Array(table::DataTable)
    return hcat(table.columns...)    
end


function Base.getindex(table::DataTable, idxs::BitArray{1})
    @assert length(idxs)==length(table.columns[1])
    subtable = DataTable(table.header)
    for i in 1:length(table.columns)
        subtable.columns[i] = table.columns[i][idxs]
    end
    return subtable
end

function Base.getindex(table::DataTable, rowindex::Int, colindex::Int)
    return table.columns[colindex][rowindex]
end

function Base.getindex(table::DataTable, rowindex::Int, colon::Colon)
    row = []
    for j=1:length(table.header)
        push!(row, table.columns[j][rowindex])
    end

    return row
end

function Base.lastindex(table::DataTable)
    length(table.columns)==0 && error("DataTable: use of 'end' in an empty table")
    return length(table.columns[1])
end


mutable struct DataBook
    tables::Array{DataTable, 1}
    function DataBook()
        this = new()
        this.tables = DataTable[]
        return this
    end
end


function push!(book::DataBook, table::DataTable)
    push!(book.tables, table)
end


function Base.getindex(book::DataBook, index::Int)
    return book.tables[index]
end


function Base.lastindex(book::DataBook)
    return length(book.tables)
end


function Base.length(book::DataBook)
    return length(book.tables)
end

function Base.iterate(book::DataBook, state=(nothing,1) )
    table, idx = state
    if idx<=length(book.tables)
        return (book.tables[idx], (book.tables[idx+1], idx+1))
    else
        return nothing
    end
end

sprintf(fmt, args...) = @eval @sprintf($fmt, $(args...))

function compress(table::DataTable, n::Int)
    nrows = length(table.columns[1])
    nrows<=n && return table[:]

    factor = (nrows-1)/(n-1)
    idxs   = [ round(Int, 1+factor*(i-1)) for i=1:n ]
    
    subtable = DataTable(table.header)
    for i=1:length(table.columns)
        subtable.columns[i] = table.columns[i][idxs]
    end

    return subtable
end

function Base.filter(table::DataTable, expr::Expr)
    fields = get_vars(expr)
    vars   = Dict{Symbol, Float64}()
    nrows  = length(table.columns[1])
    idx    = falses(nrows)
    for i=1:nrows
        for field in fields
            vars[field] = table[field][i]
        end
        idx[i] = eval_arith_expr(expr; vars...)
    end
    return table[idx]
end

function resize(table::DataTable, n::Int=0; ratio=1.0)
    np = length(table.columns[1]) # current number of points
    if n==0
        ratio > 0.0 || error("resize: ratio should be greater than zero")
        n = max(2, round(Int, np*ratio))
    end
    n > 0 || error("resize: number of rows should be greater than zero")
    np >= 4 || error("resize: Table object should contain at least 4 rows")

    ns = np - 1                # number of spacings
    nb = floor(Int64, ns/3)    # number of bezier curves
    ds = 1.0 / nb              # spacing between curves

    newtable = DataTable(table.header)
    for (j,field) in enumerate(table.header)
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

        newtable.columns[j] = V
    end

    return newtable
end

function denoise(table::DataTable, field; fraction=0.1)
    n = length(table.columns[1]) # current number of points  
    n >= 5 || error("denise: Table object should contain at least 5 rows to denoise")

    V = zeros(n)

    for i in 2:n-1
        idx1 = max(1, i-2)
        idx2 = min(n, i+2)
        

    end

    newtable = table[:]
    
end


# TODO: Improve column width for string items
function save(table::DataTable, filename::String; verbose::Bool=true, digits::Array=[])
    suitable_formats = (".dat", ".tex")

    basename, format = splitext(filename)
    format in suitable_formats || error("DataTable: cannot save in \"$format\" format. Suitable formats $suitable_formats.")

    local f::IOStream
    try
        f  = open(filename, "w")
    catch err
        warn("DataTable: File $filename could not be opened for writing.")
        return
    end

    nc = length(table.colindex)              # number of cols
    nr = nc>0 ? length(table.columns[1]) : 0 # number of rows

    if format==".dat"
        for (i,key) in enumerate(keys(table.colindex))
            @printf(f, "%12s", key)
            print(f, i!=nc ? "\t" : "\n")
        end

        # print values
        for i=1:nr
            for j=1:nc
                item = table.columns[j][i]
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

        verbose && printstyled("  file $filename written\n", color=:cyan)
    end

    if format==".tex"
        # widths calculation
        header = keys(table.colindex)
        widths = length.(header)
        types  = eltype.(table.columns)

        if length(digits)==0
            digits = repeat([4], nc)
        end
        @assert length(digits)==nc

        for (i,col) in enumerate(table.columns)
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
        for i=1:nr
            print(f, indent^level)
            for j=1:nc
                etype = types[j]
                item = table.columns[j][i]
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


function save(book::DataBook, filename::String; verbose::Bool=true)
    basename, format = splitext(filename)
    format == ".dat" || error("DataBook: cannot save in \"$format\". Suitable format is \".dat\".")

    local f::IOStream
    try
        f = open(filename, "w")
    catch err
        warn("DataBook: File $filename could not be opened for writing.")
        return
    end

    if format==".dat"

        for (k,table) in enumerate(book.tables)

            nc = length(table.colindex)              # number of cols
            nr = nc>0 ? length(table.columns[1]) : 0 # number of rows

            # print table label
            print(f, "Table (snapshot=$(k), rows=$nr)\n")

            # print header
            for (i,key) in enumerate(keys(table.colindex))
                @printf(f, "%12s", key)
                print(f, i!=nc ? "\t" : "\n")
            end

            # print values
            for i=1:nr
                for j=1:nc
                    @printf(f, "%12.5e", table.columns[j][i])
                    print(f, j!=nc ? "\t" : "\n")
                end
            end
            print(f, "\n")
        end

        verbose && printstyled("  file $filename written\n", color=:cyan)
    end
    close(f)
    return nothing

end


function DataTable(filename::String, delim='\t')
    basename, format = splitext(filename)
    format == ".dat" || error("DataTable: cannot read \"$format\". Suitable format is \".dat\".")

    if format==".dat"
        matrix, headstr = readdlm(filename, delim, header=true, use_mmap=false)
        table = DataTable(strip.(headstr), matrix)
        return table
    end
end


function DataBook(filename::String)
    delim = "\t"
    basename, format = splitext(filename)
    format == ".dat" || error("DataBook: cannot read \"$format\". Suitable format is \".dat\".")

    f      = open(filename, "r")
    book   = DataBook()
    if format==".dat"
        lines = readlines(f)
        header_expected = false
        for (i,line) in enumerate(lines)
            items = split(line)
            length(items)==0 && strip(line)=="" && continue

            if items[1]=="Table"
                header_expected = true
                continue
            end
            if header_expected # add new table
                header = [ strip(key) for key in split(line, delim) ]
                push!(book.tables, DataTable(header))
                header_expected = false
                continue
            end

            length(book.tables) == 0 && error("DataBook: Wrong file format. Use DataTable(filename) to read a table")
            row = []
            for item in items
                try
                    char1 = item[1]
                    isnumeric(char1) || char1 in ('+','-') ?  push!(row, Meta.parse(item)) : push!(row, item)
                catch err
                    @error "DataBook: Error while reading value '$item' at line $i"
                    throw(err)
                end
            end
            push!(book.tables[end], row) # udpate newest table
        end
    end

    close(f)
    return book
end

# Functions for backwards compatibility
loadtable(filename::String, delim='\t') = DataTable(filename, delim)
loadbook(filename::String) = DataBook(filename)

# TODO: Improve display. Include column datatype
function Base.show(io::IO, table::DataTable)
    println(io)
    if length(table.columns)==0
        print(io, "DataTable()")
        return
    end
    nc = length(table.colindex)     # number of fields (columns)
    nr = length(table.columns[1])   # number of rows

    if nr==0
        print(io, "DataTable()")
        return
    end

    header = keys(table.colindex)
    types  = typeof.(getindex.(table.columns,1))

    hwidths = length.(header)
    widths  = zeros(Int, length(header))
    useformat   = falses(length(header))
    shortformat = falses(length(header))

    for (i,col) in enumerate(table.columns)
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
        for (i,col) in enumerate(table.columns)
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
    for i=1:nr
        if i>half_vrows && nr-i>=half_vrows
            i==half_vrows+1 && println(" ⋮")
            continue
        end

        print(io, " │ ")
        for j=1:nc
            etype = types[j]
            item = table.columns[j][i]
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


function Base.show(io::IO, book::DataBook)
    print(io, "DataBook (tables=$(length(book.tables))):\n")
    n = length(book.tables)
    for (k,table) in enumerate(book.tables)
        # print table label
        nitems = length(table.columns[1])
        print(io, " Table (snapshot=$(k), rows=$nitems):\n")
        str = string(table)
        k<n && print(io, str, "\n")
    end
end


randtable() = DataTable(["A","B","C"], [0:10 rand().*(sin.(0:10).+(0:10)) rand().*(cos.(0:10).+(0:10)) ])
