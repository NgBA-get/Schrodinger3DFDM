using LinearAlgebra
using BenchmarkTools

# Nx = Ny = Nz = 32 # no of grid points
Nx = 32
Ny = 32
Nz = 32

t = 5 # running time
steps = 10000
Δt = t/steps
Lx = Ly = Lz = 1

x = -Nx/2+1:Nx/2 |> collect |> k -> (k*Lx/Nx)
y = -Ny/2+1:Ny/2 |> collect |> k -> (k*Ly/Ny)
z = -Nz/2+1:Nz/2 |> collect |> k -> (k*Lz/Nz)

const δx = Lx/Nx
const δy = Ly/Ny
const δz = Lz/Nz
Nlin = 0 # Non-linearity constant

V = zeros(ComplexF32, Nx, Ny, Nz)
for k1=1:Nx, k2=1:Ny, k3=1:Nz
    V[k1, k2, k3] = 32 * (x[k1]^2 + y[k2]^2 + z[k3]^2)
end

# ψ(x, y, z, t); for now ψ(x, t)

function ∇²(f)
    f₀ = f
    f₀[2:Nx-1, :, :] = (f[3:Nx, :, :] - 2f[2:Nx-1, :, :] + f[1:Nx-2, :, :]) / δx^2
    f₀[:, 2:Ny-1, :] += (f[:, 3:Ny, :] - 2f[:, 2:Ny-1, :] + f[:, 1:Ny-2, :]) / δy^2
    f₀[:, :, 2:Nz-1] += (f[:, :, 3:Nz] - 2f[:, :, 2:Nz-1] + f[:, :, 1:Nz-2]) / δz^2
    return f₀
end

function ∂ψ_∂t(ψ)
    """
    Function which outputs the time-derivative of ψ at coordinates `(x, y, z)` at time `t`.
    """
    ψ₁ = ψ
    return (-1/128 * ∇²(ψ) + 
            (V + Nlin * norm.(ψ₁)) .* ψ₁) * -1im
end
