"""
Compute the central difference approximation of the first derivative using 2-point stencil.

# Arguments
- `u1::T`: Value at the left point of the stencil.
- `u2::T`: Value at the right point of the stencil.
- `dx::T`: Grid spacing or step size.

# Returns
- `T`: Approximation of the first derivative at the midpoint of `u1` and `u2`.

# Description
The `center_2_intdiff` function calculates the central difference approximation of the first derivative of a function using a 2-point stencil. The derivative is approximated by `(u2 - u1) / dx`, where `u1` and `u2` are the function values at adjacent grid points separated by a distance `dx`.

# Example
```julia
u1 = 1.0
u2 = 2.0
dx = 0.1
derivative = center_2_intdiff(u1, u2, dx)
println(derivative)  # Output will be 10.0
"""
center_2_intdiff(u1, u2, dx) = (u2 - u1) / dx

"""
Compute the central difference approximation of the first derivative using 4-point stencil.

# Arguments
- `u1::T`: Value at the leftmost point of the stencil.
- `u2::T`: Value at the second left point of the stencil.
- `u3::T`: Value at the second right point of the stencil.
- `u4::T`: Value at the rightmost point of the stencil.
- `dx::T`: Grid spacing or step size.

# Returns
- `T`: Approximation of the first derivative at the central point of the stencil.

# Description
The `center_4_intdiff` function calculates the central difference approximation of the first derivative of a function using a 4-point stencil. The derivative is approximated by `(u1 - 27*u2 + 27*u3 - u4) / (24*dx)`, where `u1`, `u2`, `u3`, and `u4` are function values at the stencil points with spacing `dx`.

# Example
```julia
u1 = 1.0
u2 = 2.0
u3 = 3.0
u4 = 4.0
dx = 0.1
derivative = center_4_intdiff(u1, u2, u3, u4, dx)
println(derivative)  # Output will depend on input values
```
"""
center_4_intdiff(u1, u2, u3, u4, dx) = (u1 - 27*u2 + 27*u3 - u4) / (24*dx)

"""
Compute the central difference approximation of the first derivative using 6-point stencil.

# Arguments
- `u1::T`: Value at the most left point of the stencil.
- `u2::T`: Value at the second left point of the stencil.
- `u3::T`: Value at the central point of the stencil.
- `u4::T`: Value at the second right point of the stencil.
- `u5::T`: Value at the rightmost point of the stencil.
- `u6::T`: Value at the farthest right point of the stencil.
- `dx::T`: Grid spacing or step size.

# Returns
- `T`: Approximation of the first derivative at the central point of the stencil.

# Description
The `center_6_intdiff` function calculates the central difference approximation of the first derivative of a function using a 6-point stencil. The derivative is approximated by `(-9*u1 + 125*u2 - 2250*u3 + 2250*u4 - 125*u5 + 9*u6) / (1920*dx)`, where `u1` through `u6` are function values at the stencil points with spacing `dx`.

# Example
```julia
u1 = 1.0
u2 = 2.0
u3 = 3.0
u4 = 4.0
u5 = 5.0
u6 = 6.0
dx = 0.1
derivative = center_6_intdiff(u1, u2, u3, u4, u5, u6, dx)
println(derivative)  # Output will depend on input values
"""
center_6_intdiff(u1, u2, u3, u4, u5, u6,  dx) = (-9*u1 + 125*u2 - 2250*u3 + 2250*u4 - 125*u5 + 9*u6) / (1920*dx)

"""
Compute the first derivative of a 3D array using different numerical differentiation schemes along the x-direction.

# Arguments
- `u::Array{T, 3}`: 3D array representing the function values on a grid.
- `dx::T`: Grid spacing or step size in the x-direction.
- `scheme::Int`: Differentiation scheme to use. Options are:
  - `1`: 2nd-order central difference (default).
  - `2`: 4th-order central difference.
  - `3`: 6th-order central difference.
  - `4`: 6th-order compact central difference.
- `boundary::Array{Int, 1}`: Boundary condition types for the x-direction (and other directions if applicable). Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the x-direction, computed using the selected scheme.

# Description
The `int_d1x` function calculates the first derivative of a 3D array `u` with respect to the x-direction using one of the following schemes:
- **2nd-order central difference** (`scheme=1`): Uses a simple central difference approximation.
- **4th-order central difference** (`scheme=2`): Uses a more accurate central difference approximation.
- **6th-order central difference** (`scheme=3`): Uses an even more accurate central difference approximation.
- **6th-order compact central difference** (`scheme=4`): Uses a compact scheme for higher accuracy and stability.

Boundary conditions are applied based on the provided `boundary` array and `dirichlet` values, allowing for different treatment of the boundaries depending on the problem setup.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
boundary = [1, 1, 0, 0, 0, 0]  # Neumann boundary condition in x-direction
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = int_d1x(u, dx, scheme=2, boundary=boundary, dirichlet=dirichlet)
println(derivative)
````
"""
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

"""
Compute the first derivative of a 3D array using different numerical differentiation schemes along the y-direction.

# Arguments
- `u::Array{T, 3}`: 3D array representing the function values on a grid.
- `dx::T`: Grid spacing or step size in the y-direction.
- `scheme::Int`: Differentiation scheme to use. Options are:
  - `1`: 2nd-order central difference (default).
  - `2`: 4th-order central difference.
  - `3`: 6th-order central difference.
  - `4`: 6th-order compact central difference.
- `boundary::Array{Int, 1}`: Boundary condition types for the y-direction (and other directions if applicable). Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the y-direction, computed using the selected scheme.

# Description
The `int_d1y` function calculates the first derivative of a 3D array `u` with respect to the y-direction using one of the following schemes:
- **2nd-order central difference** (`scheme=1`): Uses a simple central difference approximation.
- **4th-order central difference** (`scheme=2`): Uses a more accurate central difference approximation.
- **6th-order central difference** (`scheme=3`): Uses an even more accurate central difference approximation.
- **6th-order compact central difference** (`scheme=4`): Uses a compact scheme for higher accuracy and stability.

Boundary conditions are applied based on the provided `boundary` array and `dirichlet` values, allowing for different treatment of the boundaries depending on the problem setup.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
boundary = [0, 0, 1, 1, 0, 0]  # Neumann boundary condition in y-direction
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = int_d1y(u, dx, scheme=2, boundary=boundary, dirichlet=dirichlet)
println(derivative)
````
"""
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

"""
Compute the first derivative of a 3D array using different numerical differentiation schemes along the z-direction.

# Arguments
- `u::Array{T, 3}`: 3D array representing the function values on a grid.
- `dx::T`: Grid spacing or step size in the z-direction.
- `scheme::Int`: Differentiation scheme to use. Options are:
  - `1`: 2nd-order central difference (default).
  - `2`: 4th-order central difference.
  - `3`: 6th-order central difference.
  - `4`: 6th-order compact central difference.
- `boundary::Array{Int, 1}`: Boundary condition types for the z-direction (and other directions if applicable). Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the z-direction, computed using the selected scheme.

# Description
The `int_d1z` function calculates the first derivative of a 3D array `u` with respect to the z-direction using one of the following schemes:
- **2nd-order central difference** (`scheme=1`): Uses a simple central difference approximation.
- **4th-order central difference** (`scheme=2`): Uses a more accurate central difference approximation.
- **6th-order central difference** (`scheme=3`): Uses an even more accurate central difference approximation.
- **6th-order compact central difference** (`scheme=4`): Uses a compact scheme for higher accuracy and stability.

Boundary conditions are applied based on the provided `boundary` array and `dirichlet` values, allowing for different treatment of the boundaries depending on the problem setup.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
boundary = [0, 0, 0, 0, 1, 1]  # Neumann boundary condition in z-direction
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = int_d1z(u, dx, scheme=2, boundary=boundary, dirichlet=dirichlet)
println(derivative)
````
"""
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

"""
Compute the first derivative of a 3D array using a 2nd-order central difference scheme along the x-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the x-direction.
- `b::Array{Int, 1}`: Boundary condition types for the x-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the x-direction computed using a 2nd-order central difference scheme.

# Description
The `center_int_2_d1x` function calculates the first derivative of a 3D array `u` with respect to the x-direction using a 2nd-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_2_d1x(u, dx, b, dirichlet)
println(derivative)
````
"""
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

"""
Compute the first derivative of a 3D array using a 4th-order central difference scheme along the x-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the x-direction.
- `b::Array{Int, 1}`: Boundary condition types for the x-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the x-direction computed using a 4th-order central difference scheme.

# Description
The `center_int_4_d1x` function calculates the first derivative of a 3D array `u` with respect to the x-direction using a 4th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_4_d1x(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 6th-order central difference scheme along the x-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the x-direction.
- `b::Array{Int, 1}`: Boundary condition types for the x-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the x-direction computed using a 6th-order central difference scheme.

# Description
The `center_int_6_d1x` function calculates the first derivative of a 3D array `u` with respect to the x-direction using a 6th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_6_d1x(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 2nd-order central difference scheme along the y-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the y-direction.
- `b::Array{Int, 1}`: Boundary condition types for the y-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the y-direction computed using a 2nd-order central difference scheme.

# Description
The `center_int_2_d1y` function calculates the first derivative of a 3D array `u` with respect to the y-direction using a 2nd-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_2_d1y(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 4th-order central difference scheme along the y-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the y-direction.
- `b::Array{Int, 1}`: Boundary condition types for the y-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the y-direction computed using a 4th-order central difference scheme.

# Description
The `center_int_4_d1y` function calculates the first derivative of a 3D array `u` with respect to the y-direction using a 4th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_4_d1y(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 6th-order central difference scheme along the y-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the y-direction.
- `b::Array{Int, 1}`: Boundary condition types for the y-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the y-direction computed using a 6th-order central difference scheme.

# Description
The `center_int_6_d1y` function calculates the first derivative of a 3D array `u` with respect to the y-direction using a 6th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_6_d1y(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 2nd-order central difference scheme along the z-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the z-direction.
- `b::Array{Int, 1}`: Boundary condition types for the z-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the z-direction computed using a 2nd-order central difference scheme.

# Description
The `center_int_2_d1z` function calculates the first derivative of a 3D array `u` with respect to the z-direction using a 2nd-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_2_d1z(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 4th-order central difference scheme along the z-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the z-direction.
- `b::Array{Int, 1}`: Boundary condition types for the z-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the z-direction computed using a 4th-order central difference scheme.

# Description
The `center_int_4_d1z` function calculates the first derivative of a 3D array `u` with respect to the z-direction using a 4th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_4_d1z(u, dx, b, dirichlet)
println(derivative)
```
"""
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

"""
Compute the first derivative of a 3D array using a 6th-order central difference scheme along the z-direction.

# Arguments
- `u::Array{T, 3}`: 3D array of function values on the grid.
- `dx::T`: Grid spacing or step size in the z-direction.
- `b::Array{Int, 1}`: Boundary condition types for the z-direction and other directions if applicable. Default is `[0, 0, 0, 0, 0, 0]`. Possible values are:
  - `0`: Periodic boundary condition.
  - `1`: Neumann boundary condition.
  - `2`: Dirichlet boundary condition.
- `dirichlet::Array{T, 1}`: Values for Dirichlet boundary conditions. Default is `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`.

# Returns
- `Array{T, 3}`: 3D array of the first derivative of `u` in the z-direction computed using a 6th-order central difference scheme.

# Description
The `center_int_6_d1z` function calculates the first derivative of a 3D array `u` with respect to the z-direction using a 6th-order central difference approximation. Boundary conditions are applied based on the provided `b` array and `dirichlet` values.

# Example
```julia
u = rand(10, 10, 10)
dx = 0.1
b = [0, 0, 0, 0, 0, 0]
dirichlet = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
derivative = center_int_6_d1z(u, dx, b, dirichlet)
println(derivative)
```
"""
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
