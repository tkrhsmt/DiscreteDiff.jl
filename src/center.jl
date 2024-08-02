
center1_2(u1, u2, dx) = (u2 - u1) / (2 * dx)
center1_4(u1, u2, u3, u4, dx) = (-u4 + 8 * u3 - 8 * u2 + u1) / (12 * dx)
center1_6(u1, u2, u3, u4, u5, u6, dx) = (u6 - 9 * u5 + 45 * u4 - 45 * u3 + 9 * u2 - u1) / (60 * dx)

function d1x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1x(u, dx, boundary, dirichlet)
    end

end

function d1y(u, dy, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1y(u, dy, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1y(u, dy, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1y(u, dy, boundary, dirichlet)
    end

end

function d1z(u, dz, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1z(u, dz, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1z(u, dz, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1z(u, dz, boundary, dirichlet)
    end

end

function center_2_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 1, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_2(ub[i-1, j, k], ub[i+1, j, k], dx)
            end
        end
    end

    return du

end

function center_4_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_4(ub[i-2, j, k], ub[i-1, j, k], ub[i+1, j, k], ub[i+2, j, k], dx)
            end
        end
    end

    return du

end

function center_6_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 3, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_6(ub[i-3, j, k], ub[i-2, j, k], ub[i-1, j, k], ub[i+1, j, k], ub[i+2, j, k], ub[i+3, j, k], dx)
            end
        end
    end

    return du

end

function center_2_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 1, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_2(ub[i, j-1, k], ub[i, j+1, k], dx)
            end
        end
    end

    return du

end

function center_4_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_4(ub[i, j-2, k], ub[i, j-1, k], ub[i, j+1, k], ub[i, j+2, k], dx)
            end
        end
    end

    return du

end

function center_6_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 3, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_6(ub[i, j-3, k], ub[i, j-2, k], ub[i, j-1, k], ub[i, j+1, k], ub[i, j+2, k], ub[i, j+3, k], dx)
            end
        end
    end

    return du

end

function center_2_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 1, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_2(ub[i, j, k-1], ub[i, j, k+1], dx)
            end
        end
    end

    return du

end

function center_4_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_4(ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k+1], ub[i, j, k+2], dx)
            end
        end
    end

    return du

end

function center_6_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 3, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center1_6(ub[i, j, k-3], ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k+1], ub[i, j, k+2], ub[i, j, k+3], dx)
            end
        end
    end

    return du

end
