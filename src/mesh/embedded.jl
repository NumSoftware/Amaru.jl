# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Remove line-joint-cells from mesh and set up line cells as embedded cells
function generate_embedded_cells!(mesh::Mesh)

    newcells = []
    id = 0
    for cell in mesh.elems
        if cell.shape.family==LINEJOINTCELL
            # link solid cell to line cells
            solid, line = cell.linked_elems
            line.linked_elems = [ solid ]
        else
            # save non joint1D cells
            id += 1
            cell.id = id  # update id
            push!(newcells, cell)
        end
    end

    # update mesh
    mesh.elems = newcells
    return mesh
end
