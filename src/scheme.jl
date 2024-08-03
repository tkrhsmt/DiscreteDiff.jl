"""
    CyclicThomas(a,b,c,d)
    特殊な三重対角行列を持つ連立方程式を解く

        [b1 c1 0     ...     a1
         a2 b2 c2    ...     0
    A=    :  :      :  :      :
         0  0  ... aN-1 bN-1 cN-1
         cN 0  ... 0    aN   bN]

    d= [d1, d2, ... dN]

    A * u = d
この時の`u`を解く

# Arguments
- `a` : 三重対角行列の1成分
- `b` : 三重対角行列の2成分
- `c` : 三重対角行列の3成分
- `d` : `A*u=d`における`d`
"""
function CyclicThomas(a, b, c, d)
    len = length(a)
    α = c[len]
    β = a[1]
    γ = -b[1]
    bb = zeros(len)
    u = zeros(len)
    bb[1] = b[1] - γ
    bb[len] = b[len] - α * β / γ
    for i in 2:len-1
        bb[i] = b[i]
    end
    y = Thomas(a, bb, c, d)

    u[1] = γ
    u[len] = α

    z = Thomas(a, bb, c, u)
    f = (y[1] + β * y[len] / γ) / (1.0 + z[1] + β * z[len] / γ)
    for i in 1:len
        u[i] = y[i] - f * z[i]
    end
    return u
end
# --------------------------------------------------
"""
    Thomas(a,b,c,d)
    三重対角行列を持つ連立方程式を解く

        [b1 c1 0     ...     0
        a2 b2 c2    ...     0
    A=    :  :      :  :      :
        0  0  ... aN-1 bN-1 cN-1
        0  0  ... 0    aN   bN]

    d= [d1, d2, ... dN]

    A * u = d
この時の`u`を解く

# Arguments
- `a` : 三重対角行列の1成分
- `b` : 三重対角行列の2成分
- `c` : 三重対角行列の3成分
- `d` : `A*u=d`における`d`
"""
function Thomas(a, b, c, d)
    e = length(a)
    gam = Array{Float64}(undef, e)
    x = Array{Float64}(undef, e)
    bet = b[1]
    x[1] = d[1] / bet

    for i = 2:e
        gam[i] = c[i-1] / bet
        bet = b[i] - a[i] * gam[i]
        x[i] = (d[i] - a[i] * x[i-1]) / bet
    end

    for i = e-1:-1:1
        x[i] = x[i] - gam[i+1] * x[i+1]
    end
    return x

end
