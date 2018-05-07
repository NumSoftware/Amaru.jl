# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.getindex
import Base.keys
export DTable, DBook, push!, keys, getindex, save, loadtable, loadbook


# DTable object

type DTable
    data    ::Array{Array{Float64,1},1}
    colindex::Dict{Symbol,Int} # Data index
    fields  ::Array{Symbol,1}
    function DTable()
        this = new()
        this.data     = [ ]
        this.colindex = Dict{Symbol,Int}() 
        this.fields   = []
        return this
    end
    function DTable(header::Array{Symbol,1})
        this = new()
        this.data     = [ Float64[] for s in header ]
        this.colindex = Dict( key=>i for (i,key) in enumerate(header) )
        this.fields   = copy(header)
        return this
    end
    function DTable(header::Array{Symbol,1}, data::Array{Array{Float64,1},1})
        this      = new()
        nfields   = length(header)
        ncols     = length(data)
        nfields  != ncols && error("DTable: header and data fields do not match")
        this.data     = deepcopy(data)
        this.colindex = Dict( key=>i for (i,key) in enumerate(header) )
        this.fields   = copy(header)
        return this
    end
    function DTable(header::Array{Symbol,1}, matrix::Array{Float64,2})
        this      = new()
        nfields   = length(header)
        ncols     = size(matrix,2)
        nfields  != ncols && error("DTable: header and data fields do not match")
        this.data     = [ matrix[:,i] for i=1:nfields ]
        this.colindex = Dict( key=>i for (i,key) in enumerate(header) )
        this.fields   = copy(header)
        return this
    end
end


type DBook
    tables::Array{DTable, 1}
    function DBook()
        this = new()
        this.tables = DTable[]
        return this
    end
end

import Base.push!
function push!(table::DTable, row::Array{Float64,1})
    @assert length(table.fields)==length(row)
    for (i,val) in enumerate(row)
        push!(table.data[i], val)
    end
end

function push!(book::DBook, table::DTable)
    push!(book.tables, table)
end

function keys(table::DTable)
    return table.fields
end

function Base.push!(table::DTable, dict::Associative{Symbol,Float64})
    if length(table.data)==0
        table.data     = [ Float64[v] for (k,v) in dict ]
        table.colindex = Dict( key=>i for (i,key) in enumerate(keys(dict)) )
        table.fields   = collect(keys(dict))
    else
        nrows = length(table.data[1])
        for (k,v) in dict
            # Add data
            colindex = get(table.colindex, k, 0)
            if colindex==0
                # add new column
                new_col = zeros(nrows)
                push!(new_col, v)
                push!(table.data, new_col)
                push!(table.fields, k)
                table.colindex[k] = length(table.data)
            else
                push!(table.data[colindex], v)
            end
        end
        # Add zero for missing values if any
        for col in table.data
            if length(col)==nrows
                push!(col, 0.0)
            end
        end
    end
end

function Base.getindex(table::DTable, key::Symbol)
    return table.data[table.colindex[key]]
end

function Base.getindex(table::DTable, keys::Array{Symbol,1})
    data = [ table[key] for key in keys ]
    subtable = DTable(keys, data)
    return subtable
end

function Base.getindex(book::DBook, index::Int)
    return book.tables[index]
end

function Base.endof(book::DBook)
    return length(book.tables)
end

Base.start(book::DBook) = 1
Base.next(book::DBook, idx::Int) = book.tables[idx], idx+1
Base.done(book::DBook, idx::Int) = idx == length(book.tables)

function save(table::DTable, filename::String; verbose::Bool=true)
    format = split(filename*".", ".")[2]
    format != "dat" && error("save DTable: filename should have \"dat\" extension")

    f  = open(filename, "w")
    nc = length(table.fields)     # number of fields (columns)
    nr = length(table.data[1])  # number of rows

    if format=="dat"
        # print header
        for i=1:nc
            @printf(f, "%12s", table.fields[i])
            print(f, i!=nc? "\t" : "\n")
        end

        # print values
        for i=1:nr
            for j=1:nc
                @printf(f, "%12.6e", table.data[j][i])
                print(f, j!=nc? "\t" : "\n")
            end
        end

        verbose && print_with_color(:green, "  file $filename written\n")
    end

    if format=="json"
        # enconding
        str  = JSON.json(table.colindex, 4)
        print(f, str)

        verbose && print_with_color(:green, "  file $filename written (DTable)\n")
    end

    close(f)
    return nothing
end


function save(book::DBook, filename::String; verbose::Bool=true)
    format = split(filename*".", ".")[2]
    format != "dat" && error("save DBook: filename should have \"dat\" extension")

    f  = open(filename, "w") 

    if format=="json"
        # generate dictionary
        dict_arr = [ table.colindex for table in book.tables ]
        str  = JSON.json(dict_arr, 4)
        print(f, str)

        if verbose  print_with_color(:green, "  file $filename written (DBook)\n") end
    end

    if format=="dat"

        for (k,table) in enumerate(book.tables)
            # print table label
            nitems = length(table.data[1])
            print(f, "Table (snapshot=$k, items=$nitems)\n")

            nc = length(table.colindex)     # number of fields (columns)
            nr = length(table.data[1])  # number of rows

            # print header
            for i=1:nc
                @printf(f, "%12s", table.fields[i])
                print(f, i!=nc? "\t" : "\n")
            end

            # print values
            for i=1:nr
                for j=1:nc
                    @printf(f, "%12.6e", table.data[j][i])
                    print(f, j!=nc? "\t" : "\n")
                end
            end
            print(f, "\n")
        end

        verbose && print_with_color(:green, "  file $filename written\n")
    end
    close(f)
    return nothing

end


function loadtable(filename::String)
    format = split(filename*".", ".")[2]
    format != "dat" && error("loadtable: filename should have \"dat\" extension")

    if format=="dat"
        data, headstr = readdlm(filename, '\t', header=true, use_mmap=false)
        fields = Symbol[ Symbol(strip(field)) for field in vec(headstr) ]

        table = DTable(fields , data)
        return table
    end
end

function loadbook(filename::String)
    format = split(filename*".", ".")[2]
    format != "dat" && error("loadbook: filename should have \"dat\" extension")

    f      = open(filename, "r")
    book   = DBook()
    if format=="dat"
        lines = readlines(f)
        header_expected = false
        for line in lines
            strip(line)=="" && continue
            items = split(line)
            if items[1]=="Table"
                header_expected = true
                continue
            end
            if header_expected # add new table
                fields = [ Symbol(key) for key in split(line) ]
                push!(book.tables, DTable(fields))
                header_expected = false
                continue
            end

            length(book.tables) == 0 && error("loadbook: Wrong file format. Use loadtable() to read a table")
            row = parse.(items)
            push!(book.tables[end], row) # udpate newest table
        end
    end

    close(f)
    return book
end

function Base.show(io::IO, table::DTable)
    nc = length(table.colindex)  # number of fields (columns)
    nr = length(table.data[1])   # number of rows
    matrix = [ table.data[j][i] for i=1:nr, j=1:nc ]
    header = reshape(Text.([:row; table.fields]), (1, nc+1))
    print(io, "DTable $(nr)x$nc:\n")
    Base.showarray(io, [header; collect(1:nr) matrix], false, header=false)
end

function Base.show(io::IO, book::DBook)
    print(io, "DBook (tables=$(length(book.tables))):\n")
    for (k,table) in enumerate(book.tables)
        # print table label
        nitems = length(table.data[1])
        #print(io, "  Table (snapshot=$k, items=$nitems):\n")
        str = "  "*replace(string(table), "\n","\n  ")
        print(io, str, "\n")
    end
end
