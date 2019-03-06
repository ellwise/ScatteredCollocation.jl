
function griddedSquarePoints(dx::Number)
    num_intervals = convert(Int, cld(1, dx))
    dx = 1/num_intervals
    x = [[dx*i, dx*j] for i = 0:num_intervals, j = 0:num_intervals]
    bnd = falses(num_intervals+1, num_intervals+1)
    bnd[1, :] .= true
    bnd[:, 1] .= true
    bnd[end, :] .= true
    bnd[:, end] .= true
    x = x[:] # flatten
    bnd = bnd[:]
    return x, bnd
end

function haltonSquarePoints(dx::Number)
    length = 1.0
    num_points = convert(Int, cld(length^2, dx^2))
    x = [[halton_sequence(i, 2), halton_sequence(i, 3)] for i = 1:num_points]
    bnd = [false for i = 1:num_points]
    return x, bnd
end

# below seems to basically work, but I made it up...
function haltonDiscPoints(dx::Number)
    radius = 1.0
    num_points = convert(Int, cld(pi*radius^2, dx^2))
    thetas = 2*pi*halton_sequence.(collect(1:num_points), 2)
    radii = radius*sqrt.(halton_sequence.(collect(1:num_points), 3))
    x = [radii[i] .* [sin(thetas[i]), cos(thetas[i])] for i=1:num_points]
    bnd = [false for i=1:num_points]
    return x, bnd
end

function concentricDiscPoints(dx::Number)
    PACKING_NUMBER = 2*pi
    num_radial = convert(Int, cld(1, dx))
    r_arr = dx:dx:num_radial*dx
    x = [[0.0,0.0]]
    bnd = [false]
    num_current = 0
    for i = 1:num_radial
        num_current += PACKING_NUMBER
        dtheta = 2*pi/round(Int, num_current)
        thetas = range(0, step=dtheta, length=round(Int, num_current))
        x_new = [r_arr[i].*[sin(t), cos(t)] for t=thetas]
        if i == num_radial
            bnd_new = [true for t=thetas]
        else
            bnd_new = [false for t=thetas]
        end
        x = vcat(x, x_new)
        bnd = vcat(bnd, bnd_new)
    end
    return x, bnd
end

function randomDiscPoints(dx::Number)
    radius = 1.0
    num_nodes = convert(Int, cld(pi*radius^2, dx^2))
    thetas = 2*pi*rand(num_nodes)
    radii = radius*sqrt.(rand(num_nodes))
    x = [radii[i] .* [sin(thetas[i]), cos(thetas[i])] for i=1:num_nodes]
    bnd = [false for i=1:num_nodes]
    return x, bnd
end
