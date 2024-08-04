
"""
    compact_6_d1x(u, dx, boundary, dirichlet)

Calculate the first derivative in the x direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the x direction.
"""
function compact_6_d1x(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nx)

    α = 1 / 3
    a = 14 / 9 / 2 / dx
    b = 1 / 9 / 4 / dx

    for k in 1:nz
        for j in 1:ny
            Aa = ones(nx) * α
            Ab = ones(nx)
            Ac = ones(nx) * α

            for i in 1:nx
                B[i] = a * (ub[i+1, j, k] - ub[i-1, j, k]) + b * (ub[i+2, j, k] - ub[i-2, j, k])
            end

            if bo[1] == 0 #周期境界条件
                du[:, j, k] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[1] == 2 # ディリクレ条件
                    B[1] -= center1_2(ub[-1, j, k], ub[1, j, k], dx) * α
                end
                if bo[2] == 2 # ディリクレ条件
                    B[nx] -= center1_2(ub[nx, j, k], ub[nx+2, j, k], dx) * α
                end
                du[:, j, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end


"""
    compact_6_d1y(u, dy, boundary, dirichlet)

Calculate the first derivative in the y direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the y direction.
"""
function compact_6_d1y(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(ny)

    α = 1 / 3
    a = 14 / 9 / 2 / dx
    b = 1 / 9 / 4 / dx

    for k in 1:nz
        for i in 1:nx
            Aa = ones(ny) * α
            Ab = ones(ny)
            Ac = ones(ny) * α

            for j in 1:ny
                B[j] = a * (ub[i, j+1, k] - ub[i, j-1, k]) + b * (ub[i, j+2, k] - ub[i, j-2, k])
            end

            if bo[3] == 0 #周期境界条件
                du[i, :, k] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[3] == 2 # ディリクレ条件
                    B[1] -= center1_2(ub[i, -1, k], ub[i, 1, k], dx) * α
                end
                if bo[4] == 2 # ディリクレ条件
                    B[ny] -= center1_2(ub[i, ny, k], ub[i, ny+2, k], dx) * α
                end
                du[i, :, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end


"""
    compact_6_d1z(u, dz, boundary, dirichlet)

Calculate the first derivative in the z direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the z direction.
"""
function compact_6_d1z(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nz)

    α = 1 / 3
    a = 14 / 9 / 2 / dx
    b = 1 / 9 / 4 / dx

    for j in 1:ny
        for i in 1:nx
            Aa = ones(nz) * α
            Ab = ones(nz)
            Ac = ones(nz) * α

            for k in 1:nz
                B[k] = a * (ub[i, j, k+1] - ub[i, j, k-1]) + b * (ub[i, j, k+2] - ub[i, j, k-2])
            end

            if bo[5] == 0 #周期境界条件
                du[i, j, :] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[5] == 2 # ディリクレ条件
                    B[1] -= center1_2(ub[i, j, -1], ub[i, j, 1], dx) * α
                end
                if bo[6] == 2 # ディリクレ条件
                    B[nz] -= center1_2(ub[i, j, nz], ub[i, j, nz+2], dx) * α
                end
                du[i, j, :] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end

"""
    center_2_d2x(u, dx, boundary, dirichlet)

Calculate the second derivative in the x direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the x direction.
"""
function compact_6_d2x(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nx)

    α = 2 / 11
    a = 12 / 11 / (dx^2)
    b = 3 / 11 / 4 / (dx^2)

    for k in 1:nz
        for j in 1:ny
            Aa = ones(nx) * α
            Ab = ones(nx)
            Ac = ones(nx) * α

            for i in 1:nx
                B[i] = a * (ub[i+1, j, k] - 2*ub[i, j, k] + ub[i-1, j, k]) + b * (ub[i+2, j, k] - 2*ub[i, j, k] + ub[i-2, j, k])
            end

            if bo[1] == 0 #周期境界条件
                du[:, j, k] = CyclicThomas(Aa, Ab, Ac, B)
            else #その他の境界条件
                B[1] -= center2_2(ub[-1, j, k], ub[0, j, k], ub[1, j, k], dx) * α
                B[nx] -= center2_2(ub[nx, j, k], ub[nx+1, j, k], ub[nx+2, j, k], dx) * α
                du[:, j, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end

"""
    center_2_d2y(u, dx, boundary, dirichlet)

Calculate the second derivative in the y direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the y direction.
"""
function compact_6_d2y(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(ny)

    α = 2 / 11
    a = 12 / 11 / (dx^2)
    b = 3 / 11 / 4 / (dx^2)

    for k in 1:nz
        for i in 1:nx
            Aa = ones(ny) * α
            Ab = ones(ny)
            Ac = ones(ny) * α

            for j in 1:ny
                B[j] = a * (ub[i, j+1, k] - 2*ub[i, j, k] + ub[i, j-1, k]) + b * (ub[i, j+2, k] - 2*ub[i, j, k] + ub[i, j-2, k])
            end

            if bo[3] == 0 #周期境界条件
                du[i, :, k] = CyclicThomas(Aa, Ab, Ac, B)
            else #その他の境界条件
                B[1] -= center2_2(ub[i, -1, k], ub[i, 0, k], ub[i, 1, k], dx) * α
                B[ny] -= center2_2(ub[i, ny, k], ub[i, ny+1, k], ub[i, ny+2, k], dx) * α
                du[i, :, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end

"""
    center_2_d2z(u, dx, boundary, dirichlet)

Calculate the second derivative in the z direction using 6th order compact difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the z direction.
"""
function compact_6_d2z(u, dx, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nz)

    α = 1 / 3
    a = 14 / 9 / 2 / dx
    b = 1 / 9 / 4 / dx

    for j in 1:ny
        for i in 1:nx
            Aa = ones(nz) * α
            Ab = ones(nz)
            Ac = ones(nz) * α

            for k in 1:nz
                B[k] = a * (ub[i, j, k+1] -2*ub[i, j, k] + ub[i, j, k-1]) + b * (ub[i, j, k+2] - 2*ub[i, j, k] + ub[i, j, k-2])
            end

            if bo[5] == 0 #周期境界条件
                du[i, j, :] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[5] == 2 # ディリクレ条件
                    B[1] -= center2_2(ub[i, j, -1], ub[i, j, 0], ub[i, j, 1], dx) * α
                end
                if bo[6] == 2 # ディリクレ条件
                    B[nz] -= center2_2(ub[i, j, nz], ub[i, j, nz+1], ub[i, j, nz+2], dx) * α
                end
                du[i, j, :] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end
