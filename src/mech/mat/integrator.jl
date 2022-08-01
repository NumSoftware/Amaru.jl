export MechIntegrator
export stress_update

mutable struct MechIntegrator
    ipd::IpState
    mat::Material
    table::DataTable
    function MechIntegrator(mat::Material)
        this = new()
        env = ModelEnv()
        this.ipd = ip_state_type(mat)(env)
        this.mat = mat
        this.table = DataTable()
        return this
    end
end


function stress_update(int::MechIntegrator, Δε; nincs=1)
    Δεi = Δε/nincs
    for i in 1:nincs
        stress_update(int.mat, int.ipd, Δεi)
        vals = ip_state_vals(int.mat, int.ipd)
        push!(int.table, vals)
    end
end

