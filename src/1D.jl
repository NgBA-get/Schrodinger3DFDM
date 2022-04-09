using LinearAlgebra
using Interpolations
using Plots

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
const Nlin = 0 # Non-linearity constant

V = zeros(ComplexF32, Nx, Ny, Nz)
for k1=1:Nx, k2=1:Ny, k3=1:Nz
    V[k1, k2, k3] = 32 * (x[k1]^2 + y[k2]^2 + z[k3]^2)
end

# ψ(x, y, z, t); for now ψ(x, t)

function ∇²(f, t)
    f₀ = f(:, :, :, t)
    f₀[2:Nx-1, :, :] = (f(3:Nx, :, :, t) - 2f(2:Nx-1, :, :, t) + f(1:Nx-2, :, : , t)) / δx^2
    f₀[:, 2:Ny-1, :] += (f(:, 3:Ny, :, t) - 2f(:, 2:Ny-1, :, t) + f(:, 1:Ny-2, :, t)) / δy^2
    f₀[:, :, 2:Nz-1] += (f(:, :, 3:Nz, t) - 2f(:, :, 2:Nz-1, t) + f(:, :, 1:Nz-2,t)) / δz^2
    return f₀
end

function ∂ψ_∂t(ψ, t)
    """
    Function which outputs the time-derivative of ψ at coordinates `(x, y, z)` at time `t`.
    """
    ψ₁ = ψ(:, :, :, t)
    return (-1/128 * ∇²(ψ, t) + 
            (V + Nlin * norm.(ψ₁)) .* ψ₁) * -1im
end
