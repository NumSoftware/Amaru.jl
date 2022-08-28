# Remove output files from tests

path  = dirname(@__FILE__)
#files = readdir(path)

println("  deleting temporary files...")
for (root, dirs, files) in walkdir(path)
    for file in files
        ext = split(file*".", ".")[2]
        if ext in ["vtk", "vtu", "dat", "log"]
            fullname = joinpath(root, file)
            rm(fullname)
        end
    end
end