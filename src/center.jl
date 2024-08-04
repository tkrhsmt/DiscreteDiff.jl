"""
    center1_2(u1, u2, dx)

Calculate the first derivative using 2nd order central difference.

# Arguments
- `u1::Number`: The value at the previous point.
- `u2::Number`: The value at the next point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The first derivative.
"""
center1_2(u1, u2, dx) = (u2 - u1) / (2 * dx)

"""
    center1_4(u1, u2, u3, u4, dx)

Calculate the first derivative using 4th order central difference.

# Arguments
- `u1::Number`: The value at the previous 3rd point.
- `u2::Number`: The value at the previous point.
- `u3::Number`: The value at the next point.
- `u4::Number`: The value at the next 3rd point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The first derivative.
"""
center1_4(u1, u2, u3, u4, dx) = (-u4 + 8 * u3 - 8 * u2 + u1) / (12 * dx)

"""
    center1_6(u1, u2, u3, u4, u5, u6, dx)

Calculate the first derivative using 6th order central difference.

# Arguments
- `u1::Number`: The value at the previous 5th point.
- `u2::Number`: The value at the previous 3rd point.
- `u3::Number`: The value at the previous point.
- `u4::Number`: The value at the next point.
- `u5::Number`: The value at the next 3rd point.
- `u6::Number`: The value at the next 5th point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The first derivative.
"""
center1_6(u1, u2, u3, u4, u5, u6, dx) = (u6 - 9 * u5 + 45 * u4 - 45 * u3 + 9 * u2 - u1) / (60 * dx)

"""
    center2_2(u1, u2, u3, dx)

Calculate the second derivative using 2nd order central difference.

# Arguments
- `u1::Number`: The value at the previous point.
- `u2::Number`: The value at the current point.
- `u3::Number`: The value at the next point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The second derivative.
"""
center2_2(u1, u2, u3, dx) = (u3 + -2*u2 + u1) / (dx^2)

"""
    center2_4(u1, u2, u3, u4, u5, dx)

Calculate the second derivative using 4th order central difference.

# Arguments
- `u1::Number`: The value at the previous 2nd point.
- `u2::Number`: The value at the previous point.
- `u3::Number`: The value at the current point.
- `u4::Number`: The value at the next point.
- `u5::Number`: The value at the next 2nd point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The second derivative.
"""
center2_4(u1, u2, u3, u4, u5, dx) = (- u5 + 16 * u4 - 30 * u3 + 16 * u2 - u1) / (12 * dx^2)

"""
    center2_6(u1, u2, u3, u4, u5, u6, u7, dx)

Calculate the second derivative using 6th order central difference.

# Arguments
- `u1::Number`: The value at the previous 3rd point.
- `u2::Number`: The value at the previous 2nd point.
- `u3::Number`: The value at the previous point.
- `u4::Number`: The value at the current point.
- `u5::Number`: The value at the next point.
- `u6::Number`: The value at the next 2nd point.
- `u7::Number`: The value at the next 3rd point.
- `dx::Number`: The distance between points.

# Returns
- `Number`: The second derivative.
"""
center2_6(u1, u2, u3, u4, u5, u6, u7, dx) = (2 * u7 - 27 * u6 + 270 * u5 - 490 * u4 + 270 * u3 - 27 * u2 + 2 * u1) / (180 * dx^2)

"""
    d1x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the first derivative in the x direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order, 4: 6th order compact).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the x direction.
"""
function d1x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1x(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d1x(u, dx, boundary, dirichlet)
    end

end

"""
    d1y(u, dy, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the first derivative in the y direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order, 4: 6th order compact).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the y direction.
"""
function d1y(u, dy, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1y(u, dy, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1y(u, dy, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1y(u, dy, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d1y(u, dy, boundary, dirichlet)
    end

end

"""
    d1z(u, dz, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the first derivative in the z direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order, 4: 6th order compact).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the z direction.
"""
function d1z(u, dz, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d1z(u, dz, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d1z(u, dz, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d1z(u, dz, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d1z(u, dz, boundary, dirichlet)
    end

end

"""
    d2x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the second derivative in the x direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the x direction.
"""
function d2x(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d2x(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d2x(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d2x(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d2x(u, dx, boundary, dirichlet)
    end

end

"""
    d2y(u, dy, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the second derivative in the y direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the y direction.
"""
function d2y(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d2y(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d2y(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d2y(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d2y(u, dx, boundary, dirichlet)
    end

end

"""
    d2z(u, dz, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

Calculate the second derivative in the z direction using specified scheme.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `scheme::Int`: The scheme to use (1: 2nd order, 2: 4th order, 3: 6th order).
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the z direction.
"""
function d2z(u, dx, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心差分
        return center_2_d2z(u, dx, boundary, dirichlet)
    elseif scheme == 2 #4次中心差分
        return center_4_d2z(u, dx, boundary, dirichlet)
    elseif scheme == 3 #6次中心差分
        return center_6_d2z(u, dx, boundary, dirichlet)
    elseif scheme == 4 #6次compact差分
        return compact_6_d2z(u, dx, boundary, dirichlet)
    end

end

"""
    center_2_d1x(u, dx, boundary, dirichlet)

Calculate the first derivative in the x direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the x direction.
"""
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

"""
    center_4_d1x(u, dx, boundary, dirichlet)

Calculate the first derivative in the x direction using 4th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the x direction.
"""
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


"""
    center_6_d1x(u, dx, boundary, dirichlet)

Calculate the first derivative in the x direction using 6th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the x direction.
"""
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

"""
    center_2_d1y(u, dy, boundary, dirichlet)

Calculate the first derivative in the y direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the y direction.
"""
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

"""
    center_4_d1y(u, dy, boundary, dirichlet)

Calculate the first derivative in the y direction using 4th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the y direction.
"""
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

"""
    center_6_d1y(u, dy, boundary, dirichlet)

Calculate the first derivative in the y direction using 6th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dy::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the y direction.
"""
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

"""
    center_2_d1z(u, dz, boundary, dirichlet)

Calculate the first derivative in the z direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the z direction.
"""
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

"""
    center_4_d1z(u, dz, boundary, dirichlet)

Calculate the first derivative in the z direction using 4th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the z direction.
"""
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

"""
    center_6_d1z(u, dz, boundary, dirichlet)

Calculate the first derivative in the z direction using 6th order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dz::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The first derivative in the z direction.
"""
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

"""
    center_2_d2x(u, dx, boundary, dirichlet)

Calculate the second derivative in the x direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the x direction.
"""
function center_2_d2x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 1, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_2(ub[i-1, j, k], ub[i, j, k], ub[i+1, j, k], dx)
            end
        end
    end

    return du

end

"""
    center_4_d2x(u, dx, boundary, dirichlet)

Calculate the second derivative in the x direction using 4nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the x direction.
"""
function center_4_d2x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_4(ub[i-2, j, k], ub[i-1, j, k], ub[i, j, k], ub[i+1, j, k], ub[i+2, j, k], dx)
            end
        end
    end

    return du

end

"""
    center_6_d2x(u, dx, boundary, dirichlet)

Calculate the second derivative in the x direction using 6nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the x direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the x direction.
"""
function center_6_d2x(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 3, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_6(ub[i-3, j, k], ub[i-2, j, k], ub[i-1, j, k], ub[i, j, k], ub[i+1, j, k], ub[i+2, j, k], ub[i+3, j, k], dx)
            end
        end
    end

    return du

end

"""
    center_2_d2y(u, dx, boundary, dirichlet)

Calculate the second derivative in the y direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the y direction.
"""
function center_2_d2y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 1, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_2(ub[i, j-1, k], ub[i, j, k], ub[i, j+1, k], dx)
            end
        end
    end

    return du

end

"""
    center_4_d2y(u, dx, boundary, dirichlet)

Calculate the second derivative in the y direction using 4nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the y direction.
"""
function center_4_d2y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_4(ub[i, j-2, k], ub[i, j-1, k], ub[i, j, k], ub[i, j+1, k], ub[i, j+2, k], dx)
            end
        end
    end

    return du

end

"""
    center_6_d2y(u, dx, boundary, dirichlet)

Calculate the second derivative in the y direction using 6nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the y direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the y direction.
"""
function center_6_d2y(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 3, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_6(ub[i, j-3, k], ub[i, j-2, k], ub[i, j-1, k], ub[i, j, k], ub[i, j+1, k], ub[i, j+2, k], ub[i, j+3, k], dx)
            end
        end
    end

    return du

end

"""
    center_2_d2z(u, dx, boundary, dirichlet)

Calculate the second derivative in the z direction using 2nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the z direction.
"""
function center_2_d2z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 1, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_2(ub[i, j, k-1], ub[i, j, k], ub[i, j, k+1], dx)
            end
        end
    end

    return du

end

"""
    center_4_d2z(u, dx, boundary, dirichlet)

Calculate the second derivative in the x direction using 4nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the z direction.
"""
function center_4_d2z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_4(ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k], ub[i, j, k+1], ub[i, j, k+2], dx)
            end
        end
    end

    return du

end

"""
    center_6_d2z(u, dx, boundary, dirichlet)

Calculate the second derivative in the z direction using 6nd order central difference.

# Arguments
- `u::Array`: The input 3D array.
- `dx::Number`: The distance between points in the z direction.
- `boundary::Array`: The boundary conditions.
- `dirichlet::Array`: The values for Dirichlet boundary conditions.

# Returns
- `Array`: The second derivative in the z direction.
"""
function center_6_d2z(u, dx, b, dirichlet)

    nx, ny, nz = size(u)
    du = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 3, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                du[i, j, k] = center2_6(ub[i, j, k-3], ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k], ub[i, j, k+1], ub[i, j, k+2], ub[i, j, k+3], dx)
            end
        end
    end

    return du

end
