center_2_intdiff(u1, u2, dx) = (u2 - u1) / dx
center_4_intdiff(u1, u2, u3, u4, dx) = (u1 - 27*u2 + 27*u3 - u4) / (24*dx)
center_6_intdiff(u1, u2, u3, u4, u5, u6,  dx) = (-9*u1 + 125*u2 - 2250*u3 + 2250*u4 - 125*u5 + 9*u6) / (1920*dx)


function int_d1x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間差分
        return center_int_2_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間差分
        return center_int_4_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間差分
        return center_int_6_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間差分
        return compact_6_int_d1x(u, dx, boundary, dirichlet)
    end
end

function int_d1y(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間差分
        return center_int_2_d1y(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間差分
        return center_int_4_d1y(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間差分
        return center_int_6_d1y(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間差分
        return compact_6_int_d1y(u, dx, boundary, dirichlet)
    end
end

function int_d1z(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間差分
        return center_int_2_d1z(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間差分
        return center_int_4_d1z(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間差分
        return center_int_6_d1z(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間差分
        return compact_6_int_d1z(u, dx, boundary, dirichlet)
    end
end

function center_int_2_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 1, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_2_intdiff(ub[i, j, k], ub[i+1, j, k], dx)
            end
        end
    end

    return intdu
end

function center_int_4_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_4_intdiff(ub[i-1, j, k], ub[i, j, k],ub[i+1, j, k], ub[i+2, j, k], dx)
            end
        end
    end

    return intdu
end

function center_int_6_d1x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 3, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_6_intdiff(ub[i-2, j, k], ub[i-1, j, k],ub[i, j, k], ub[i+1, j, k],ub[i+2, j, k], ub[i+3, j, k], dx)
            end
        end
    end

    return intdu
end

function center_int_2_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 1, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_2_intdiff(ub[i, j, k], ub[i, j+1, k], dx)
            end
        end
    end

    return intdu
end

function center_int_4_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_4_intdiff(ub[i, j-1, k], ub[i, j, k],ub[i, j+1, k], ub[i, j+2, k], dx)
            end
        end
    end

    return intdu
end

function center_int_6_d1y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 3, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_6_intdiff(ub[i, j-2, k], ub[i, j-1, k],ub[i, j, k], ub[i, j+1, k],ub[i, j+2, k], ub[i, j+3, k], dx)
            end
        end
    end

    return intdu
end

function center_int_2_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 1, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_2_intdiff(ub[i, j, k], ub[i, j, k+1], dx)
            end
        end
    end

    return intdu
end

function center_int_4_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_4_intdiff(ub[i, j, k-1], ub[i, j, k],ub[i, j, k+1], ub[i, j, k+2], dx)
            end
        end
    end

    return intdu
end

function center_int_6_d1z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    intdu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 3, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intdu[i, j, k] = center_6_intdiff(ub[i, j, k-2], ub[i, j, k-1],ub[i, j, k], ub[i, j, k+1],ub[i, j, k+2], ub[i, j, k+3], dx)
            end
        end
    end

    return intdu
end
