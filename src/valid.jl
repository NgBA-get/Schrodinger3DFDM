using Trapz

function I_nsqre(grid)
    """
    grid: 3D array of ψ
    """
    return trapz((x, y, z), grid .|> k->norm(k)^2)
end
