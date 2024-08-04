"""
    center_2_int(u1, u2)

Calculate the integral using 2nd order central difference.

# Arguments
- `u1::Number`: The value at the first point.
- `u2::Number`: The value at the second point.

# Returns
- `Number`: The integral value.
"""
center_2_int(u1, u2) = (u1 + u2) / 2.0

"""
    center_4_int(u1, u2, u3, u4)

Calculate the integral using 4th order central difference.

# Arguments
- `u1::Number`: The value at the first point.
- `u2::Number`: The value at the second point.
- `u3::Number`: The value at the third point.
- `u4::Number`: The value at the fourth point.

# Returns
- `Number`: The integral value.
"""
center_4_int(u1, u2, u3, u4) = (-u1 + 9*u2 + 9*u3 - u4) / 16.0

"""
    center_6_int(u1, u2, u3, u4, u5, u6)

Calculate the integral using 6th order central difference.

# Arguments
- `u1::Number`: The value at the first point.
- `u2::Number`: The value at the second point.
- `u3::Number`: The value at the third point.
- `u4::Number`: The value at the fourth point.
- `u5::Number`: The value at the fifth point.
- `u6::Number`: The value at the sixth point.

# Returns
- `Number`: The integral value.
"""
center_6_int(u1, u2, u3, u4, u5, u6) = (3*u1 - 25*u2 + 150*u3 + 150*u4 - 25*u5 + 3*u6) / 256.0

"""
Compute interpolated values using different schemes.

# Arguments
- `u::Array{T}`: Input array of values to be interpolated. Type `T` should be compatible with the interpolation schemes.
- `scheme::Int`: Interpolation scheme to use. The options are:
  - `1`: 2nd-order central interpolation
  - `2`: 4th-order central interpolation
  - `3`: 6th-order central interpolation
  - `4`: 6th-order compact interpolation
- `boundary::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain. Default is `[0, 0, 0, 0, 0, 0]`.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- The interpolated values as per the selected scheme.

# Description
The `intx` function applies interpolation to the input data `u` using one of the specified schemes. The choice of scheme determines the order and type of interpolation used:
- Scheme 1 uses 2nd-order central interpolation.
- Scheme 2 uses 4th-order central interpolation.
- Scheme 3 uses 6th-order central interpolation.
- Scheme 4 uses 6th-order compact interpolation.

The function also takes into account boundary conditions specified by the `boundary` and `dirichlet` arrays.
"""
function intx(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_x(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_x(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_x(u, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間
        return compact_6_intx(u, boundary, dirichlet)
    end
end

"""
Compute interpolated values along the y-direction using different schemes.

# Arguments
- `u::Array{T}`: Input array of values to be interpolated along the y-direction. Type `T` should be compatible with the interpolation schemes.
- `scheme::Int`: Interpolation scheme to use. The options are:
  - `1`: 2nd-order central interpolation
  - `2`: 4th-order central interpolation
  - `3`: 6th-order central interpolation
  - `4`: 6th-order compact interpolation
- `boundary::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain. Default is `[0, 0, 0, 0, 0, 0]`.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- The interpolated values along the y-direction as per the selected scheme.

# Description
The `inty` function applies interpolation to the input data `u` along the y-direction using one of the specified schemes. The choice of scheme determines the order and type of interpolation used:
- Scheme 1 uses 2nd-order central interpolation.
- Scheme 2 uses 4th-order central interpolation.
- Scheme 3 uses 6th-order central interpolation.
- Scheme 4 uses 6th-order compact interpolation.

The function also considers boundary conditions specified by the `boundary` and `dirichlet` arrays.
"""
function inty(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_y(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_y(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_y(u, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間
        return compact_6_inty(u, boundary, dirichlet)
    end
end

"""
Compute interpolated values along the z-direction using different schemes.

# Arguments
- `u::Array{T}`: Input array of values to be interpolated along the z-direction. Type `T` should be compatible with the interpolation schemes.
- `scheme::Int`: Interpolation scheme to use. The options are:
  - `1`: 2nd-order central interpolation
  - `2`: 4th-order central interpolation
  - `3`: 6th-order central interpolation
  - `4`: 6th-order compact interpolation
- `boundary::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain. Default is `[0, 0, 0, 0, 0, 0]`.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- The interpolated values along the z-direction as per the selected scheme.

# Description
The `intz` function applies interpolation to the input data `u` along the z-direction using one of the specified schemes. The choice of scheme determines the order and type of interpolation used:
- Scheme 1 uses 2nd-order central interpolation.
- Scheme 2 uses 4th-order central interpolation.
- Scheme 3 uses 6th-order central interpolation.
- Scheme 4 uses 6th-order compact interpolation.

The function also considers boundary conditions specified by the `boundary` and `dirichlet` arrays.
"""
function intz(u, scheme=1, boundary=[0, 0, 0, 0, 0, 0], dirichlet=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    if scheme == 1 #2次中心補間
        return int_2_z(u, boundary, dirichlet)
    elseif scheme == 2 #4次中心補間
        return int_4_z(u, boundary, dirichlet)
    elseif scheme == 3 #6次中心補間
        return int_6_z(u, boundary, dirichlet)
    elseif scheme == 4 #6次compact補間
        return compact_6_intz(u, boundary, dirichlet)
    end
end

"""
Perform 2nd-order central interpolation in the x-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the x-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the x-direction using 2nd-order central interpolation.

# Description
The `int_2_x` function applies 2nd-order central interpolation along the x-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the x-direction using the `center_2_int` function.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_2_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 1, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i, j, k], ub[i+1, j, k])
            end
        end
    end

    return intu
end

"""
Perform 2nd-order central interpolation in the y-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the y-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the y-direction using 2nd-order central interpolation.

# Description
The `int_2_y` function applies 2nd-order central interpolation along the y-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the y-direction using the `center_2_int` function.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_2_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 1, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i, j, k], ub[i, j+1, k])
            end
        end
    end

    return intu
end

"""
Perform 2nd-order central interpolation in the z-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the z-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the z-direction using 2nd-order central interpolation.

# Description
The `int_2_z` function applies 2nd-order central interpolation along the z-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the z-direction using the `center_2_int` function.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_2_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 1, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_2_int(ub[i, j, k], ub[i, j, k+1])
            end
        end
    end

    return intu
end

"""
Perform 4th-order central interpolation in the x-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the x-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the x-direction using 4th-order central interpolation.

# Description
The `int_4_x` function applies 4th-order central interpolation along the x-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the x-direction using the `center_4_int` function, which uses values from adjacent points in the x-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_4_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 2, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i-1, j, k], ub[i, j, k], ub[i+1, j, k], ub[i+2, j, k])
            end
        end
    end

    return intu
end

"""
Perform 4th-order central interpolation in the y-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the y-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the y-direction using 4th-order central interpolation.

# Description
The `int_4_y` function applies 4th-order central interpolation along the y-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the y-direction using the `center_4_int` function, which uses values from adjacent points in the y-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_4_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 2, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i, j-1, k], ub[i, j, k],ub[i, j+1, k], ub[i, j+2, k])
            end
        end
    end

    return intu
end

"""
Perform 4th-order central interpolation in the z-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the z-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the z-direction using 4th-order central interpolation.

# Description
The `int_4_z` function applies 4th-order central interpolation along the z-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the z-direction using the `center_4_int` function, which uses values from adjacent points in the z-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_4_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 2, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_4_int(ub[i, j, k-1], ub[i, j, k], ub[i, j, k+1], ub[i, j, k+2])
            end
        end
    end

    return intu
end

"""
Perform 6th-order central interpolation in the x-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the x-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the x-direction using 6th-order central interpolation.

# Description
The `int_6_x` function applies 6th-order central interpolation along the x-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the x-direction using the `center_6_int` function, which uses values from multiple adjacent points in the x-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_6_x(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 3, 0, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i-2, j, k], ub[i-1, j, k], ub[i, j, k], ub[i+1, j, k], ub[i+2, j, k], ub[i+3, j, k])
            end
        end
    end

    return intu
end

"""
Perform 6th-order central interpolation in the y-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the y-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the y-direction using 6th-order central interpolation.

# Description
The `int_6_y` function applies 6th-order central interpolation along the y-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the y-direction using the `center_6_int` function, which uses values from multiple adjacent points in the y-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_6_y(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 3, 0, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i, j-2, k], ub[i, j-1, k],ub[i, j, k], ub[i, j+1, k],ub[i, j+2, k], ub[i, j+3, k])
            end
        end
    end

    return intu
end

"""
Perform 6th-order central interpolation in the z-direction.

# Arguments
- `u::Array{T}`: Input 3D array of values to be interpolated. Type `T` should be compatible with the interpolation function.
- `b::Array{T, 1}`: Array of boundary conditions. Length should be 6, corresponding to the boundary conditions on each edge of the domain in the z-direction.
- `dirichlet::Array{T, 1}`: Array of Dirichlet boundary conditions. Length should be 6.

# Returns
- `intu::Array{T}`: Interpolated 3D array along the z-direction using 6th-order central interpolation.

# Description
The `int_6_z` function applies 6th-order central interpolation along the z-direction to the input data `u`. The function first adds boundary conditions to the input array `u` using the `add_boundary` function. It then performs central interpolation at each point in the z-direction using the `center_6_int` function, which uses values from multiple adjacent points in the z-direction for higher accuracy.

Boundary conditions specified by `b` and Dirichlet boundary conditions specified by `dirichlet` are applied to handle the edges of the domain.
"""
function int_6_z(u, b, dirichlet)

    nx, ny, nz = size(u)
    intu = zeros(nx, ny, nz)

    ub = add_boundary(u, 0, 0, 3, b[1], b[2], b[3], b[4], b[5], b[6], dirichlet)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                intu[i, j, k] = center_6_int(ub[i, j, k-2], ub[i, j, k-1], ub[i, j, k], ub[i, j, k+1], ub[i, j, k+2], ub[i, j, k+3])
            end
        end
    end

    return intu
end
