

"""
    compact_6_d1x(u,dx)
    差分法を利用したx方向1階微分

境界条件：周期境界
差分法：6次compact差分

# Arguments
- `u` : 微分を求める配列
- `dx` : x方向の分割幅
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
    compact_6_d1y(u,dx)
    差分法を利用したy方向1階微分

境界条件：周期境界
差分法：6次compact差分

# Arguments
- `u` : 微分を求める配列
- `dx` : y方向の分割幅
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
    compact_6_d1z(u,dx)
    差分法を利用したz方向1階微分

境界条件：周期境界
差分法：6次compact差分

# Arguments
- `u` : 微分を求める配列
- `dx` : z方向の分割幅
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
