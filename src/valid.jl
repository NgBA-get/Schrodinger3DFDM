using Trapz

function grid_norm(grid)
    return trapz((x, y, z), grid .|> k->norm(k)^2)
end

# function energy(grid)
    
# end
