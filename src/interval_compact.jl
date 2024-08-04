
function compact_6_intx(u, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 4, 0, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nx)

    α = 0.461658
    d = 0.00146508
    a = (75 + 70 * α - 640 * d) / 128
    b = (126 * α - 25 + 2304 * d) / 256
    c = (3 - 10 * α - 1280 * d) / 256

    for k in 1:nz
        for j in 1:ny
            Aa = ones(nx) * α
            Ab = ones(nx)
            Ac = ones(nx) * α

            for i in 1:nx
                B[i] = (a * (ub[i+1, j, k] + ub[i, j, k])
                        + b * (ub[i+2, j, k] + ub[i-1, j, k])
                        + c * (ub[i+3, j, k] + ub[i-2, j, k])
                        + d * (ub[i+4, j, k] + ub[i-3, j, k]))
            end

            if bo[1] == 0 #周期境界条件
                du[:, j, k] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[1] == 1 #ノイマン条件
                    B[1] -= ub[0, j, k] * α
                elseif bo[1] == 2 # ディリクレ条件
                    B[1] -= dirichlet[1] * α
                end
                if bo[2] == 1 #ノイマン条件
                    B[1] -= ub[nx+1, j, k] * α
                elseif bo[2] == 2 # ディリクレ条件
                    B[nx] -= dirichlet[2] * α
                end
                du[:, j, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end


function compact_6_inty(u, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 4, 0, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(ny)

    α = 0.461658
    d = 0.00146508
    a = (75 + 70 * α - 640 * d) / 128
    b = (126 * α - 25 + 2304 * d) / 256
    c = (3 - 10 * α - 1280 * d) / 256

    for k in 1:nz
        for i in 1:nx
            Aa = ones(ny) * α
            Ab = ones(ny)
            Ac = ones(ny) * α

            for j in 1:ny
                B[j] = (a * (ub[i, j+1, k] + ub[i, j, k])
                        + b * (ub[i, j+2, k] + ub[i, j-1, k])
                        + c * (ub[i, j+3, k] + ub[i, j-2, k])
                        + d * (ub[i, j+4, k] + ub[i, j-3, k]))
            end

            if bo[3] == 0 #周期境界条件
                du[i, :, k] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[3] == 1 #ノイマン条件
                    B[1] -= ub[i, 0, k] * α
                elseif bo[3] == 2 # ディリクレ条件
                    B[1] -= dirichlet[3] * α
                end
                if bo[4] == 1 #ノイマン条件
                    B[1] -= ub[i, ny+1, k] * α
                elseif bo[4] == 2 # ディリクレ条件
                    B[nx] -= dirichlet[4] * α
                end
                du[i, :, k] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end


function compact_6_intz(u, bo, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 4, bo[1], bo[2], bo[3], bo[4], bo[5], bo[6], dirichlet)

    B = zeros(nz)

    α = 0.461658
    d = 0.00146508
    a = (75 + 70 * α - 640 * d) / 128
    b = (126 * α - 25 + 2304 * d) / 256
    c = (3 - 10 * α - 1280 * d) / 256

    for j in 1:ny
        for i in 1:nx
            Aa = ones(nz) * α
            Ab = ones(nz)
            Ac = ones(nz) * α

            for k in 1:nz
                B[k] = (a * (ub[i, j, k+1] + ub[i, j, k])
                        + b * (ub[i, j, k+2] + ub[i, j, k-1])
                        + c * (ub[i, j, k+3] + ub[i, j, k-2])
                        + d * (ub[i, j, k+4] + ub[i, j, k-3]))
            end

            if bo[5] == 0 #周期境界条件
                du[i, j, :] = CyclicThomas(Aa, Ab, Ac, B)
            else
                if bo[5] == 1 #ノイマン条件
                    B[1] -= ub[i, j, 0] * α
                elseif bo[5] == 2 # ディリクレ条件
                    B[1] -= dirichlet[5] * α
                end
                if bo[6] == 1 #ノイマン条件
                    B[1] -= ub[i, j, nz+1] * α
                elseif bo[6] == 2 # ディリクレ条件
                    B[nx] -= dirichlet[6] * α
                end
                du[i, j, :] = Thomas(Aa, Ab, Ac, B)
            end

        end
    end

    return du
end
