# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

# Remove joint cells from mesh and set up line cells as embedded cells
function generate_embedded_cells!(mesh::Mesh)

    newcells = []
    id = 0
    for cell in mesh.cells
        if cell.shape.family==JOINT1D_SHAPE
            # link solid cell to line cells
            solid, line = cell.linked_cells
            line.linked_cells = [ solid ]
        else
            # save non joint1D cells
            id += 1
            cell.id = id  # update id
            push!(newcells, cell)
        end
    end

    # update mesh
    mesh.cells = newcells
    return mesh
end
