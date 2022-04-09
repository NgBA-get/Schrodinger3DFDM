include("1D.jl")

ϕ = Array{Complex{Float32}, 4}(undef, Nx, Ny, Nz, steps+1) # Keeping operations in Float32 to keep array size low and speed high.

# Initializing ψ(x, y, z, 1) at t = 1 (as array indexing starts from 1 in Julia).
# Initializing to a Gaussian Normal Distribution.

gauss_3d(x, y, z) = (64/π)^(3/4)*exp(-64/2*(x^2 + y^2 + z^2)) + 0im |> Complex{Float32}
gauss_2d(x, y) = (64/π)^(1/2)*exp(-64/2(x^2 + y^2)) + 0im |> Complex{Float32}

for idx in CartesianIndices(ϕ[:, :, :, 1])
    xx, yy, zz = idx.I
    ϕ[idx.I..., 1] = gauss_3d(x[xx], y[yy], z[zz])
end
# 3D gaussian initialization done.

ψ(x, y, z, t) = ϕ[x, y, z, t] # To keep up consistency with the rest of the codebase.

function ψ_set(val, t)
    ϕ[:, :, :, t] = val
    return nothing
end
