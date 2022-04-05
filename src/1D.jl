using LinearAlgebra
using Interpolations

Nx = Ny = Nz = 32 # no of grid points
t = 10 # running time
steps = 10000
Δt = (t-1)/steps
Lx = Ly = Lz = 1

x = -Nx/2+1:Nx/2 |> collect |> k -> (k*Lx/Nx)
y = -Ny/2+1:Ny/2 |> collect |> k -> (k*Ly/Ny)
z = -Nz/2+1:Nz/2 |> collect |> k -> (k*Lz/Nz)

const δx = Lx/Nx
const δy = Ly/Ny
const δz = Lz/Nz
const Nlin = 0 # Non-linearity constant

# (x, y, z) are grid points, not actual euclidean space coordinates.

function V(x...)
    """
    Function representing velocity potential.
    Calculates the potential for at a point in Rⁿ space.
    """
    return 64/2*mapreduce(k->k^2, +, x)
end

# ψ(x, y, z, t); for now ψ(x, t)

function ∇²(f, x, y, z, t)
    """
    Calculates the Laplacian of the function `f = f(x, y, z, t)`. 
    
    Arguments:
        f: function of which we need the Laplacian
        x: grid point in x-direction
        y: grid point in y-direction
        z: grid point in z-direction
        t: time coordinate
    """
    dx = (f(x+1, y, z, t) - 2f(x, y, z, t) + f(x-1, y, z ,t)) / δx^2
    dy = (f(x, y+1, z, t) - 2f(x, y, z, t) + f(x, y-1, z ,t)) / δy^2
    dz = (f(x, y, z+1, t) - 2f(x, y, z, t) + f(x, y, z-1 ,t)) / δz^2
    
    return dx+dy+dz
end

function ∂ψ_∂t(ψ, V, x, y, z, t)
    """
    Function which outputs the time-derivative of ψ at coordinates `(x, y, z)` at time `t`.
    """
    return (-1/128 * ∇²(ψ, x, y, z, t) + 
            V(x, y, z, t) * ψ(x, y, z, t) +
            Nlin * norm(ψ(x, y, z, t)) * ψ(x, y, z, t)) * -1im
end
