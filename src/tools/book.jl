# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DataBook, push!, save, loadbook


# DataBook type
mutable struct DataBook
    tables::Array{DataTable, 1}
    function DataBook()
        this = new()
        this.tables = DataTable[]
        return this
    end
end


Base.length(book::DataBook) = length(book.tables)
Base.size(book::DataBook)   = (length(book.tables),)


function Base.push!(book::DataBook, table::DataTable)
    push!(book.tables, table)
end


function Base.getindex(book::DataBook, index::Int)
    return book.tables[index]
end


function Base.lastindex(book::DataBook)
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

            columns = getcolumns(table)
            colidx  = getcolidx(table)
            header  = getheader(table)
            nr, nc  = size(table)

            # nc = length(colidx)              # number of cols
            # nr = nc>0 ? length(columns[1]) : 0 # number of rows

            # print table label
            print(f, "Table (snapshot=$(k), rows=$nr)\n")

            # print header
            for (i,key) in enumerate(keys(colidx))
                @printf(f, "%12s", key)
                print(f, i!=nc ? "\t" : "\n")
            end

            # print values
            for i=1:nr
                for j=1:nc
                    @printf(f, "%12.5e", columns[j][i])
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


# function for backwards compatibility
loadbook(filename::String) = DataBook(filename)


function Base.show(io::IO, book::DataBook)
    print(io, "DataBook (tables=$(length(book.tables))):\n")
    n = length(book.tables)
    for (k,table) in enumerate(book.tables)
        # print table label
        nitems = length(columns[1])
        print(io, " Table (snapshot=$(k), rows=$nitems):\n")
        str = string(table)
        k<n && print(io, str, "\n")
    end
end