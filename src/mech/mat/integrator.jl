export MechIntegrator
export update_state

mutable struct MechIntegrator
    state::IpState
    mat::Material
    table::DataTable
    function MechIntegrator(mat::Material)
        this = new()
        ctx = Context()
        this.state = compat_state_type(mat)(ctx)
        this.mat = mat
        this.table = DataTable()
        return this
    end
end


function update_state!(int::MechIntegrator, Δε; nincs=1)
    Δεi = Δε/nincs
    for i in 1:nincs
        update_state!(int.mat, int.state, Δεi)
        vals = ip_state_vals(int.mat, int.state)
        push!(int.table, vals)
    end
end

