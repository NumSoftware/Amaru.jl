# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

import Base.getindex
import Base.keys
export DTable, DBook, push!, keys, getindex, save, loadtable, loadbook


# DTable object
const KeyType = Union{Symbol,AbstractString}
#const ItemType  = Union{Int64,Float64,AbstractString}
const ColType = Array{T,1} where T

mutable struct DTable
    columns ::Array{ColType,1}
    colindex::OrderedDict{String,Int} # Data index
    header::Array{String,1}
    function DTable()
        this = new()
        this.columns  = []
        this.colindex = OrderedDict() 
        this.header   = []
        return this
    end
    function DTable(header::Array)
        this = new()
        header = vec(header)
        this.columns  = [ [] for s in header ]
        this.colindex = OrderedDict( string(key)=>i for (i,key) in enumerate(header) )
        this.header   = string.(header)
        return this
    end
end


function DTable(header::Array, columns::Array{<:ColType,1})
    this      = DTable(header)
    nfields   = length(header)
    ncols     = length(columns)
    nfields  != ncols && error("DTable: header and number of data columns do not match")
    this.columns = deepcopy(columns)
    return this
end


function DTable(header::Array, matrix::Array{T,2} where T)
    this   = DTable(header)
    nkeys  = length(header)
    ncols  = size(matrix,2)
    nkeys != ncols && error("DTable: header and number of data columns do not match")
    types = [ typeof(matrix[1,i]) for i=1:ncols ]
    this.columns = [ convert(Array{types[i],1}, matrix[:,i]) for i=1:ncols ]
    return this
end


mutable struct DBook
    tables::Array{DTable, 1}
    function DBook()
        this = new()
        this.tables = DTable[]
        return this
    end
end


import Base.push!
function push!(table::DTable, row::Array{T,1} where T)
    @assert length(table.colindex)==length(row)
    
    if length(table.columns[1])==0
        table.columns = [ typeof(v)[v] for v in row  ]
    else
        for (i,val) in enumerate(row)
            @show typeof(table.columns[i])
            push!(table.columns[i], val)
        end
    end
end


function push!(book::DBook, table::DTable)
    push!(book.tables, table)
end

function keys(table::DTable)
    return keys(table.colindex)
end

function Base.push!(table::DTable, dict::AbstractDict)
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

function Base.getindex(table::DTable, key::KeyType)
    return table.columns[table.colindex[string(key)]]
end

function Base.getindex(table::DTable, keys::Array{<:KeyType,1})
    columns = [ table[string(key)] for key in keys ]
    subtable = DTable(keys, columns)
    return subtable
end

function Base.getindex(table::DTable, rowindex::Int, colon::Colon)
    row = []
    for j=1:length(table.header)
        push!(row, table.columns[j][rowindex])
    end

    return row
end

function Base.lastindex(table::DTable, idx::Int)
    length(table.columns)==0 && error("DTable: use of 'end' in an empty table")
    return length(table.columns[1])
end

function Base.getindex(book::DBook, index::Int)
    return book.tables[index]
end

function Base.lastindex(book::DBook)
    return length(book.tables)
end

# TODO: Check this function
function Base.iterate(book::DBook, state=(nothing,1) )
    table, idx = state
    if idx<=length(book.tables)
        return (book.tables[idx], (book.tables[i+1], idx+1))
    else
        return nothing
    end
end

sprintf(fmt, args...) = @eval @sprintf($fmt, $(args...))

# TODO: Improve column width for string items
#function save(table::DTable, filename::String; verbose::Bool=true, digits::Array{Int,1}=[])
function save(table::DTable, filename::String; verbose::Bool=true, digits::Array=[])
    format = split(filename, ".")[end]
    format in ("dat","json","tex") || error("save DTable: invalid extension $format")

    local f::IOStream
    try
        f  = open(filename, "w") 
    catch err
        @warn "DTable: File $filename could not be opened for writing."
        return
    end


    nc = length(table.colindex)             # number of cols
    nr = nc>0 ? length(table.columns[1]) : 0 # number of rows

    if format=="dat"
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

    if format=="json"
        # enconding
        str  = JSON.json(table.colindex, 4)
        print(f, str)

        verbose && printstyled("  file $filename written (DTable)\n", color=:cyan)
    end
    
    if format=="tex"
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
        println(f, raw"\begin{tabular}{", "c"^nc, "}" )
        println(f, raw"    \toprule")
        print(f, "    ")
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
        println(f, raw"    \hline")
        for i=1:nr
            print(f, "    ")
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
                    print(f, rpad(item, width))
                else
                    str = string(item)
                    print(f, rpad(item, width))
                end
                j<nc && print(f, " & ")
            end
            println(f, raw" \\\\")
        end
        println(f, raw"    \bottomrule")

        # printing ending
        println(f, raw"\end{tabular}")
    end

    close(f)
    return nothing
end


function save(book::DBook, filename::String; verbose::Bool=true)
    format = split(filename, ".")[end]
    format != "dat" && error("save DBook: filename should have \"dat\" extension")

    local f::IOStream
    try
        f  = open(filename, "w") 
    catch err
        @warn "DBook: File $filename could not be opened for writing."
        return
    end

    if format=="json"
        # generate dictionary
        dict_arr = [ table.colindex for table in book.tables ]
        str  = JSON.json(dict_arr, 4)
        print(f, str)

        if verbose  printstyled("  file $filename written (DBook)\n", color=:cyan) end
    end

    if format=="dat"

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


function loadtable(filename::String, delim='\t')
    format = split(filename, ".")[end]
    format != "dat" && error("loadtable: filename should have \"dat\" extension")

    if format=="dat"
        matrix, headstr = readdlm(filename, delim, header=true, use_mmap=false)
        table = DTable(strip.(headstr), matrix)
        return table
    end
end


function loadbook(filename::String)
    delim = "\t"
    format = split(filename, ".")[end]
    format != "dat" && error("loadbook: filename should have \"dat\" extension")

    f      = open(filename, "r")
    book   = DBook()
    if format=="dat"
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
                header = [ key for key in split(line, delim) ]
                push!(book.tables, DTable(header))
                header_expected = false
                continue
            end

            length(book.tables) == 0 && error("loadbook: Wrong file format. Use loadtable() to read a table")
            row = []
            for item in items
                try
                    char1 = item[1]
                    isnumeric(char1) || char1 in ('+','-') ?  push!(row, Meta.parse(item)) : push!(row, item)
                    #isnumeric(item[1]) ?  push!(row, Meta.parse(item)) : push!(row, item)
                catch err
                    @error "loadbook: Error while reading value '$item' at line $i"
                    throw(err)
                end
            end
            push!(book.tables[end], row) # udpate newest table
        end
    end

    close(f)
    return book
end

# TODO: Improve display. Include column datatype
function Base.show(io::IO, table::DTable)
    if length(table.columns)==0
        print("DTable()")
        return
    end
    nc = length(table.colindex)     # number of fields (columns)
    nr = length(table.columns[1])   # number of rows

    if nr==0
        print("DTable()")
        return
    end

    header = keys(table.colindex)
    widths = length.(header) 
    types  = eltype.(table.columns)
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

    print(" │ ")
    for (i,key) in enumerate(header)
        etype = types[i]
        width = widths[i]
        if etype<:Real
            print(lpad(key, width))
        else
            print(rpad(key, width))
        end
        print(" │ ")
    end
    println()

    visible_rows = 30
    half_vrows = div(visible_rows,2)

    # print values
    for i=1:nr
        if i>half_vrows && nr-i>=half_vrows
            i==half_vrows+1 && println(" ⋮")
            continue
        end

        print(" │ ")
        for j=1:nc
            etype = types[j]
            item = table.columns[j][i]
            if etype<:AbstractFloat
                item = @sprintf("%12.5e", item)
                print(lpad(item, widths[j]))
            elseif etype<:Integer
                item = @sprintf("%6d", item) 
                print(lpad(item, widths[j]))
            elseif etype<:AbstractString
                str = item
                print(rpad(item, widths[j]))
            else
                str = string(item)
                print(rpad(item, widths[j]))
            end
            print(" │ ")
        end
        i<nr && println()
    end

end

function Base.show(io::IO, book::DBook)
    print(io, "DBook (tables=$(length(book.tables))):\n")
    n = length(book.tables)
    for (k,table) in enumerate(book.tables)
        # print table label
        nitems = length(table.columns[1])
        print(io, " Table (snapshot=$(k), rows=$nitems):\n")
        str = string(table)
        k<n && print(io, str, "\n")
    end
end
