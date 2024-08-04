center_2_int(u1, u2) = (u1 + u2) / 2.0
center_4_int(u1, u2, u3, u4) = (-u1 + 9*u2 + 9*u3 - u4) / 16.0
center_6_int(u1, u2, u3, u4, u5, u6) = (-u1 + 9*u2 - 45*u3 + 45*u4 - 9*u5 + u6) / 64.0

function intx(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_x(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_x(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_x(u, boundary, dirichlet)
    end
end

function inty(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_y(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_y(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_y(u, boundary, dirichlet)
    end
end

function intz(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_z(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_z(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_z(u, boundary, dirichlet)
    end
end

function int_2_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 1, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i-1, j, k], ub[i+1, j, k])
            end
        end
    end

    return intu
end

function int_2_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 1, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i, j-1, k], ub[i, j+1, k])
            end
        end
    end

    return intu
end

function int_2_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 1, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i, j, k-1], ub[i, j, k+1])
            end
        end
    end

    return intu
end

function int_4_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i-2, j, k], ub[i-1, j, k], ub[i+1, j, k], ub[i+2, j, k])
            end
        end
    end

    return intu
end

function int_4_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i, j-2, k], ub[i, j-1, k],ub[i, j+1, k], ub[i, j+2, k])
            end
        end
    end

    return intu
end

function int_4_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k+1], ub[i, j, k+2])
            end
        end
    end

    return intu
end

function int_6_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 3, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i-3, j, k], ub[i-2, j, k], ub[i-1, j, k], ub[i+1, j, k], ub[i+2, j, k], ub[i+3, j, k])
            end
        end
    end

    return intu
end

function int_6_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 3, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i, j-3, k], ub[i, j-2, k],ub[i, j-1, k], ub[i, j+1, k],ub[i, j+2, k], ub[i, j+3, k])
            end
        end
    end

    return intu
end

function int_6_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 3, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i, j, k-3], ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k+1], ub[i, j, k+2], ub[i, j, k+3])
            end
        end
    end

    return intu
end
