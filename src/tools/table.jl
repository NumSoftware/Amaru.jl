# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.getindex
import Base.keys
export DTable, DBook, push!, keys, getindex, save, loadtable, loadbook


# DTable object

type DTable
    data  ::Array{Array{Float64,1},1}
    dict  ::Dict{Symbol,Array{Float64,1}} # Data index
    function DTable()
        this = new()
        this.dict = Dict{Symbol,Array{Float64,1}}() 
        return this
    end
    function DTable(header::Array{Symbol,1}, matrix::Array{Float64,2}=zeros(0,0))
        this = new()
        if length(matrix)==0
            this.data = [ [] for s in header]
        else
            nh = length(header)
            nf = size(matrix,2)
            if nh != nf; error("DTable: header and data fields do not match") end
            this.data = [ matrix[:,i] for i=1:nh]
        end
        this.dict = Dict( k=>v for (k,v) in zip(header, this.data) )
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
function push!(table::DTable, row::Array{Float64})
    for (i,val) in enumerate(row)
        push!(table.data[i], val)
    end
end

function push!(book::DBook, table::DTable)
    push!(book.tables, table)
end

function keys(table::DTable)
    return keys(table.dict)
end

function push!(table::DTable, dict::Dict{Symbol,Float64})
    if length(table.dict)==0
        table.data   = [ [v] for (k,v) in dict ]
        table.dict   = Dict( k=>v for (k,v) in zip(keys(dict), table.data) )
    else
        nrows = length(table.data[1])
        for (k,v) in dict
            # Add data
            if haskey(table.dict, k)
                push!(table[k], v)
            else
                # add new column
                new_arr = zeros(nrows)
                push!(new_arr, v)
                push!(table.data, new_arr)
                table.dict[k] = new_arr
            end
        end
        # Add zero for missing values if any
        for arr in table.data
            if length(arr)==nrows
                push!(arr, 0.0)
            end
        end
    end
end

function getindex(table::DTable, field::Symbol)
    return table.dict[field]
end

function getindex(book::DBook, index::Int64)
    return book.tables[index]
end


function save(table::DTable, filename::String; verbose::Bool=true)
    format = split(filename*".", ".")[2]
    f  = open(filename, "w")
    nc = length(table.dict)     # number of fields (columns)
    nr = length(table.data[1])  # number of rows

    if format=="dat"
        # map for ordered header
        ord_header = sort(collect(keys(table.dict)))
        ord_table  = [ table[key] for key in ord_header ]

        # print header
        for i=1:nc
            @printf(f, "%18s", ord_header[i])
            print(f, i!=nc? "\t" : "\n")
        end

        # print values
        for i=1:nr
            for j=1:nc
                @printf(f, "%18.10e", ord_table[j][i])
                print(f, j!=nc? "\t" : "\n")
            end
        end

        verbose && print_with_color(:green, "  file $filename written\n")
    end

    if format=="json"
        # enconding
        str  = JSON.json(table.dict, 4)
        print(f, str)

        verbose && print_with_color(:green, "  file $filename written (DTable)\n")
    end

    close(f)
    return nothing
end


function save(book::DBook, filename::String; verbose::Bool=true)
    format = split(filename*".", ".")[2]
    f  = open(filename, "w") 

    if format=="json"
        # generate dictionary
        dict_arr = [ table.dict for table in book.tables ]
        str  = JSON.json(dict_arr, 4)
        print(f, str)

        if verbose  print_with_color(:green, "  file $filename written (DBook)\n") end
    end

    if format=="dat"

        for (k,table) in enumerate(book.tables)
            # print table label
            nitems = length(table.data[1])
            print(f, "Table (snapshot=$k, items=$nitems)\n")

            nc = length(table.dict)     # number of fields (columns)
            nr = length(table.data[1])  # number of rows

            ord_header = sort(collect(keys(table.dict)))
            ord_table  = [ table[key] for key in ord_header ]

            # print header
            for i=1:nc
                @printf(f, "%18s", ord_header[i])
                print(f, i!=nc? "\t" : "\n")
            end

            # print values
            for i=1:nr
                for j=1:nc
                    @printf(f, "%18.10e", ord_table[j][i])
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
    if format=="dat"
        data, headstr = readdlm(filename, '\t', header=true, use_mmap=false)
        header = Symbol[ Symbol(strip(field)) for field in vec(headstr) ]

        table = DTable(header, data)
        return table
    end
end

function loadbook(filename::String)
    format = split(filename*".", ".")[2]
    f      = open(filename, "r")
    book   = DBook()
    if format=="dat"
        while !eof(f)
            line = readline(f)
            line=="" && continue
            items = split(line)
            if items[1]=="Table"
                keys = [ Symbol(key) for key in split(readline(f)) ]
                push!(book.tables, DTable(keys))
            else
                length(book.tables) == 0 && error("loadbook: Wrong file format. Use loadtable() to read a table")
                row = parse.(items)
                push!(book.tables[end], row)
            end
        end
    end

    close(f)
    return book
end
