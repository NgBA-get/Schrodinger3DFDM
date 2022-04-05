include("drive.jl")

function solve(ψ, cdts, t)
    """
    cdts::CartesianIndices : set of all coordinates we want to solve for.
    Solves for time instance `t`.
    """
    for cd in cdts
        k₁ = Δt * ∂ψ_∂t(ψ, V, cd.I..., t)
        k₂ = Δt * (∂ψ_∂t(ψ, V, cd.I..., t) + k₁/2)
        ψ_set(ψ(cd.I..., t) + k₂, cd.I..., t+1)
    end
end

function solveAll(ψ, cdts)
    tst = 1:steps |> collect
    for ts in tst
        solve(ψ, cdts, ts)
        println("Solved for timestep=$(ts * Δt)")
    end
    return nothing
end

# Values are exploding
