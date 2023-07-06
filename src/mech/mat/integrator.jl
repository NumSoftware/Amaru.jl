export MechIntegrator
export update_state

mutable struct MechIntegrator
    state::IpState
    matparams::MatParams
    table::DataTable
    function MechIntegrator(matparams::MatParams)
        this = new()
        env = ModelEnv()
        this.state = ip_state_type(matparams)(env)
        this.matparams = matparams
        this.table = DataTable()
        return this
    end
end


function update_state(int::MechIntegrator, Δε; nincs=1)
    Δεi = Δε/nincs
    for i in 1:nincs
        update_state(int.matparams, int.state, Δεi)
        vals = ip_state_vals(int.matparams, int.state)
        push!(int.table, vals)
    end
end

