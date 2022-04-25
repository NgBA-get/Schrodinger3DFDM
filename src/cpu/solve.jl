include("drive.jl")

function solve()
    ψ₀ = ϕ
    ψ₁ = ∂ψ_∂t(ϕ)
    ϕ_set(ψ₀ + ψ₁*Δt/2)
    ψ₁ = ∂ψ_∂t(ϕ)
    ϕ_set(ψ₀ + ψ₁*Δt)
end

function solveAll()
    for ts in 1:steps
        solve()
        if ts % 100 == 0
            println("Solved for timestep=$(ts * Δt)")
        end
    end
    return nothing
end
