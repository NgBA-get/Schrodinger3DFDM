using Trapz
include("drive.jl")

function solve(ψ, t)
    ψ₀ = ∂ψ_∂t(ψ, t)
    ψ_set(ψ₀*(1 + Δt + (Δt^2)/2), t+1)
end

function solveAll(ψ)
    for ts in 1:steps
        solve(ψ, ts)
        println("Solved for timestep=$(ts * Δt)")
    end
    return nothing
end
