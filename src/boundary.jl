# ==================================================
#                   Boundary.jl
#   Date    2024/7/31
#   Writer  T.Hashimoto
# ==================================================

using OffsetArrays

function add_boundary(u, adnx, adny, adnz, bx1, bxe, by1, bye, bz1, bze, dirichlet)

    #元の配列の大きさ
    nx, ny, nz = size(u)

    #境界を代入する配列の作成
    temp = zeros(nx + 2 * adnx, ny + 2 * adny, nz + 2 * adnz)
    unew = OffsetArray(temp, 1-adnx:nx+adnx, 1-adny:ny+adny, 1-adnz:nz+adnz)
    temp = nothing

    #ずらした配列の最初と最後の番号
    xb = 1 - adnx
    yb = 1 - adny
    zb = 1 - adnz
    xe = nx + adnx
    ye = ny + adny
    ze = nz + adnz

    #元の配列を代入
    unew[1:nx, 1:ny, 1:nz] = u

    # 周期境界 ----------------------------------------

    if bx1 == 0 #x方向の境界1
        # xb:0の配列に代入
        unew[xb:0, 1:ny, 1:nz] = u[nx-adnx+1:nx, 1:ny, 1:nz]
        unew[xb:0, yb:0, 1:nz] = u[nx-adnx+1:nx, ny-adny+1:ny, 1:nz]
        unew[xb:0, 1:ny, zb:0] = u[nx-adnx+1:nx, 1:ny, nz-adnz+1:nz]
        unew[xb:0, ny+1:ye, 1:nz] = u[nx-adnx+1:nx, 1:adny, 1:nz]
        unew[xb:0, 1:ny, nz+1:ze] = u[nx-adnx+1:nx, 1:ny, 1:adnz]
        unew[xb:0, yb:0, zb:0] = u[nx-adnx+1:nx, ny-adny+1:ny, nz-adnz+1:nz]
        unew[xb:0, yb:0, nz+1:ze] = u[nx-adnx+1:nx, ny-adny+1:ny, 1:adnz]
        unew[xb:0, ny+1:ye, zb:0] = u[nx-adnx+1:nx, 1:adny, nz-adnz+1:nz]
        unew[xb:0, ny+1:ye, nz+1:ze] = u[nx-adnx+1:nx, 1:adny, 1:adnz]
    end

    if bxe == 0 #x方向の境界2
        # nx+1:xeの配列に代入
        unew[nx+1:xe, 1:ny, 1:nz] = u[1:adnx, 1:ny, 1:nz]
        unew[nx+1:xe, yb:0, 1:nz] = u[1:adnx, ny-adny+1:ny, 1:nz]
        unew[nx+1:xe, 1:ny, zb:0] = u[1:adnx, 1:ny, nz-adnz+1:nz]
        unew[nx+1:xe, ny+1:ye, 1:nz] = u[1:adnx, 1:adny, 1:nz]
        unew[nx+1:xe, 1:ny, nz+1:ze] = u[1:adnx, 1:ny, 1:adnz]
        unew[nx+1:xe, yb:0, zb:0] = u[1:adnx, ny-adny+1:ny, nz-adnz+1:nz]
        unew[nx+1:xe, yb:0, nz+1:ze] = u[1:adnx, ny-adny+1:ny, 1:adnz]
        unew[nx+1:xe, ny+1:ye, zb:0] = u[1:adnx, 1:adny, nz-adnz+1:nz]
        unew[nx+1:xe, ny+1:ye, nz+1:ze] = u[1:adnx, 1:adny, 1:adnz]

    end

    if by1 == 0 #y方向の境界1
        # yb:0の配列に代入
        unew[1:nx, yb:0, 1:nz] = u[1:nx, ny-adny+1:ny, 1:nz]
        unew[1:nx, yb:0, zb:0] = u[1:nx, ny-adny+1:ny, nz-adnz+1:nz]
        unew[1:nx, yb:0, nz+1:ze] = u[1:nx, ny-adny+1:ny, 1:adnz]
    end

    if bye == 0 #y方向の境界2
        # ny+1:yeの配列に代入
        unew[1:nx, ny+1:ye, 1:nz] = u[1:nx, 1:adny, 1:nz]
        unew[1:nx, ny+1:ye, zb:0] = u[1:nx, 1:adny, nz-adnz+1:nz]
        unew[1:nx, ny+1:ye, nz+1:ze] = u[1:nx, 1:adny, 1:adnz]
    end

    if bz1 == 0 #z方向の境界1
        # zb:0の配列に代入
        unew[1:nx, 1:ny, zb:0] = u[1:nx, 1:ny, nz-adnz+1:nz]
    end

    if bze == 0 #z方向の境界2
        # zb:0の配列に代入
        unew[1:nx, 1:ny, nz+1:ze] = u[1:nx, 1:ny, 1:adnz]
    end

    # ノイマン境界 ----------------------------------------

    if bx1 == 1
        # xb:0の配列に代入
        for i in xb:0
            unew[i, 1:ny, 1:nz] = u[1, 1:ny, 1:nz]
            for j in yb:0
                unew[i, j, 1:nz] = u[1, 1, 1:nz]
                for k in zb:0
                    unew[i, j, k] = u[1, 1, 1]
                end
                for k in nz+1:ze
                    unew[i, j, k] = u[1, 1, nz]
                end
            end
            for k in zb:0
                unew[i, 1:ny, k] = u[1, 1:ny, 1]
            end
            for j in ny+1:ye
                unew[i, j, 1:nz] = u[1, ny, 1:nz]
                for k in zb:0
                    unew[i, j, k] = u[1, ny, 1]
                end
                for k in nz+1:ze
                    unew[i, j, k] = u[1, ny, nz]
                end
            end
            for k in nz+1:ze
                unew[i, 1:ny, k] = u[1, 1:ny, nz]
            end
        end

    end

    if bxe == 1
        # nx+1:xeの配列に代入
        for i in nx+1:xe
            unew[i, 1:ny, 1:nz] = u[nx, 1:ny, 1:nz]
            for j in yb:0
                unew[i, j, 1:nz] = u[nx, 1, 1:nz]
                for k in zb:0
                    unew[i, j, k] = u[nx, 1, 1]
                end
                for k in nz+1:ze
                    unew[i, j, k] = u[nx, 1, nz]
                end
            end
            for k in zb:0
                unew[i, 1:ny, k] = u[nx, 1:ny, 1]
            end
            for j in ny+1:ye
                unew[i, j, 1:nz] = u[nx, ny, 1:nz]
                for k in zb:0
                    unew[i, j, k] = u[nx, ny, 1]
                end
                for k in nz+1:ze
                    unew[i, j, k] = u[nx, ny, nz]
                end
            end
            for k in nz+1:ze
                unew[i, 1:ny, k] = u[nx, 1:ny, nz]
            end
        end

    end

    if by1 == 1
        # yb:0の配列に代入
        for j in yb:0
            unew[1:nx, j, 1:nz] = u[1:nx, 1, 1:nz]
            for k in zb:0
                unew[1:nx, j, k] = u[1:nx, 1, 1]
            end
            for k in nz+1:ze
                unew[1:nx, j, k] = u[1:nx, 1, nz]
            end
        end

    end

    if bye == 1
        # ny+1:yeの配列に代入
        for j in ny+1:ye
            unew[1:nx, j, 1:nz] = u[1:nx, ny, 1:nz]
            for k in zb:0
                unew[1:nx, j, k] = u[1:nx, ny, 1]
            end
            for k in nz+1:ze
                unew[1:nx, j, k] = u[1:nx, ny, nz]
            end
        end

    end

    if bz1 == 1
        # zb:0の配列に代入
        for k in zb:0
            unew[1:nx, 1:ny, k] = u[1:nx, 1:ny, 1]
        end
    end

    if bze == 1
        # nz+1:zeの配列に代入
        for k in nz+1:ze
            unew[1:nx, 1:ny, k] = u[1:nx, 1:ny, nz]
        end
    end

    # ディリクレ境界 ----------------------------------------

    if bx1 == 2
        # xb:0の配列に代入
        unew[xb:0, :, :] .= dirichlet[1]
    end

    if bxe == 2
        # nx+1:xeの配列に代入
        unew[nx+1:xe, :, :] .= dirichlet[2]
    end

    if by1 == 2
        # yb:0の配列に代入
        unew[:, yb:0, :] .= dirichlet[3]
    end

    if bye == 2
        # ny+1:yeの配列に代入
        unew[:, ny+1:ye, :] .= dirichlet[4]
    end

    if bz1 == 2
        # zb:0の配列に代入
        unew[:, :, zb:0] .= dirichlet[5]
    end

    if bze == 2
        # nz+1:zeの配列に代入
        unew[:, :, nz+1:ze] .= dirichlet[6]
    end

    return unew

end
