import GeneralizedSylvesterSolver: GeneralizedSylvesterWs, solvi_real_eliminate!,
    QuasiUpperTriangular, solve1!, transformation1, transform2,
    generalized_sylvester_solver!, solvi_complex_eliminate!, solviip_complex_eliminate!,
    solviip, solviip2, solviip_real_eliminate!
    
    
using Random
using Test

    
Random.seed!(123)

if false
n = 3
x = randn(n,n*n)
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
d = a*x + b*x*kron(c,c)
dold = d
ws = EyePlusAtKronBWs(n,n,n,2)

index = 1
depth = 1
nd = n^depth
t = QuasiUpperTriangular(lu(a\b)[2])
s = QuasiUpperTriangular(lu(c)[2])
d = zeros(n^2)
r = 1.0
d[1:n] = randn(n)
y = copy(d[1:n])
solvi_real_eliminate!(index,n,nd,1:nd,depth-1,r,t,s,d,ws)
d_target = -kron(s[[1],:]',t)*y
@test d[n+1:n^2] ≈ d_target[n+1:n^2]

index = 1
depth = 3
ws = EyePlusAtKronBWs(n,n,n,depth)
nd = n^depth
t = QuasiUpperTriangular(lu(a\b)[2])
s = QuasiUpperTriangular(lu(c)[2])
d = zeros(n^(depth+1))
r = 1.0
d[1:nd] = randn(nd)
y = copy(d[1:nd])
solvi_real_eliminate!(index,n,nd,1:nd,depth-1,r,t,s,d,ws)
d_target = -kron(s[1,:],kron(kron(s',s'),t))*y
@test d[nd+1:n^(depth+1)] ≈ d_target[nd+1:n^(depth+1)]

index = 2
d = zeros(n^(depth+1))
r = 1.0
d[nd+(1:nd)] = randn(nd)
y = copy(d[nd+(1:nd)])
solvi_real_eliminate!(index,n,nd,nd+(1:nd),depth-1,r,t,s,d,ws)
d_target = -kron(s[2,:],kron(s',s'),t)*y
@test d[2*nd+1:n^(depth+1)] ≈ d_target[2*nd+1:n^(depth+1)]

aa = a\b
t = QuasiUpperTriangular(lu(aa)[2])
s = QuasiUpperTriangular(lu(c)[2])

r = 1.0
d = randn(n)
d_orig = copy(d)
td = similar(d)
order = 0
depth = order
t2 = t*t
s2 = s*s
solve1!(r,depth,t,t2,s,s2,d,ws)
d_target = (I(n) + t)\d_orig
@test d_target ≈ d

r = 1.0
d = randn(n*n*n)
d_orig = copy(d)
td = similar(d)
order = 2
depth = order
solve1!(r,depth,t,t2,s,s2,d,ws)

d_target = (I(n^3) + kron(kron(s',s'),t))\d_orig

@test d ≈ d_target

a = 0.5
b1 = 0.1
b2 = 1.3
depth = 1
tt = QuasiUpperTriangular(t)
ss = QuasiUpperTriangular(s)
d = randn(2*n^(depth+1))
d_orig = copy(d)
transformation1(a,b1,b2,depth,tt,ss,d,ws)
d_target = d_orig + kron([a -b1; -b2 a],kron(s',t))*d_orig

@test d ≈ d_target

a = 0.5
b1 = 0.1
b2 = 1.3
depth = 2
tt = QuasiUpperTriangular(t)
ss = QuasiUpperTriangular(s)
d = randn(2*n^(depth+1))
d_orig = copy(d)
transformation1(a,b1,b2,depth,tt,ss,d,ws)
d_target = d_orig + kron([a -b1; -b2 a],kron(kron(s',s'),t))*d_orig

@test d ≈ d_target

index = 1
depth = 3
nd = n^depth
t = QuasiUpperTriangular(lu(aa)[2])
s = QuasiUpperTriangular(lu(c)[2])
s[1,2] = 0.5
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
d = zeros(n^(depth+1))
d[1:nd] = randn(nd)
d_orig = copy(d)
r1 = 1.0
r2 = 0.8
td = []
solviip_real_eliminate!(index,n,nd,1:nd,depth-1,r1,r2,t,t2,s,s2,d,ws)
d_target = -2*r1*kron(s[index,:],kron(kron(s',s'),t))*d_orig[1:nd] - (r1*r1+r2*r2)*kron(s2[index,:],kron(kron(s2',s2'),t2))*d_orig[1:nd]
@test d[nd+1:n^(depth+1)] ≈ d_target[nd+1:n^(depth+1)]
           
depth = 2
nd = n^depth
gamma = 0.3
delta1 = -0.4
delta2 = 0.6
d = randn(2*nd)
d_orig = copy(d)
transform2(r1, r2, gamma, delta1, delta2, nd, depth, t, t2, s, s2, d, ws)
G = [gamma delta1; delta2 gamma]
d_target = (I(2*nd) + 2*r1*kron(G,kron(s',t)) + (r1*r1 + r2*r2)*kron(G*G,kron(s2',t2)))*d_orig
@test d ≈ d_target

index = 1
depth = 3
nd = n^depth
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
d = zeros(n^(depth+1))
d[1:2*nd] = randn(2*nd)
d_orig = copy(d)
r1 = 1.0
r2 = 0.8
td = []
solviip_complex_eliminate!(index,n,nd,1:nd,depth-1,r1,r2,t,t2,s,s2,d,ws)
d_target = (-2*r1*kron(kron(s[1,:],kron(s',s')),t)*d_orig[1:nd] - (r1*r1+r2*r2)*kron(kron(s2[1,:],kron(s2',s2')),t2)*d_orig[1:nd]
            -2*r1*kron(kron(s[2,:],kron(s',s')),t)*d_orig[nd+1:2*nd] - (r1*r1+r2*r2)*kron(kron(s2[2,:],kron(s2',s2')),t2)*d_orig[nd+1:2*nd])
@test d[2*nd+1:n^3] ≈ d_target[2*nd+1:n^3]

depth = 1
n = 2
nd = n^depth
t = QuasiUpperTriangular(schur([1 2; -5 4])[1])
s = QuasiUpperTriangular(schur([1 -3; 2 3])[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
d = randn(n^(depth+1))
d_orig = copy(d)
alpha = s[1,1]
# s is transposed
beta1 = s[2,1]
beta2 = s[1,2]
beta = sqrt(-beta1*beta2)
td = similar(d)
solviip(alpha,beta,depth,t,t2,s,s2,d,ws)
d_target = (I(n^2) + 2*alpha*kron(s',t) + (alpha*alpha + beta*beta)*kron(s2',t2))\d_orig 
@test d ≈ d_target

depth = 1
n = 3
srand(1) #first diag block 2x2
#srand(2) # upper triangular
#srand(3) # second diag block 2x2
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
ws = EyePlusAtKronBWs(n,n,n,2)
nd = n^depth
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
d = randn(n^(depth+1))
d_orig = copy(d)
alpha = s[1,1]
# s is transposed
beta1 = s[2,1]
beta2 = s[1,2]
r1 = alpha
r2 = sqrt(-beta1*beta2)
td = similar(d)
solviip(r1,r2,depth,t,t2,s,s2,d,ws)
d_target = (I(n^2) + 2*alpha*kron(s',t) + (alpha*alpha + r2*r2)*kron(s2',t2))\d_orig 
@test d ≈ d_target

depth = 2
nd = n^depth
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
d = randn(n^3)
d_orig = copy(d)
alpha = s[1,1]
# s is transposed
beta1 = s[2,1]
beta2 = s[1,2]
r1 = alpha
r2 = sqrt(-beta1*beta2)
td = similar(d)
solviip(r1,r2,depth,t,t2,s,s2,d,ws)
d_target = (I(n^3) + 2*alpha*kron(s',kron(s',t)) + (alpha*alpha + r2*r2)*kron(kron(s2',s2'),t2))\d_orig
@test d ≈ d_target


n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
order = 1
d = randn(n^(order+1))
d_orig = copy(d)
ws = EyePlusAtKronBWs(n,n,n,order)
println("begin test")
solver!(t,s,d,order,ws)
d_target = (I(n^(order+1)) + kron(s',t))\d_orig
#d_target = (I(n^(order+1)) + kron(kron(s',s'),t))\d_orig
@test d ≈ d_target

n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
order = 3
d = randn(n^(order+1))
d_orig = copy(d)
ws = EyePlusAtKronBWs(n,n,n,order)
println("begin test")
solver!(t,s,d,order,ws)
d_target = (I(n^(order+1)) + kron(s',kron(kron(s',s'),t)))\d_orig
@test d ≈ d_target

n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b)[1])
s = QuasiUpperTriangular(schur(c)[1])
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
order = 3
ws = EyePlusAtKronBWs(n,n,n,order)
d = randn(n^(order+1))
d_orig = copy(d)
solver!(t,s,d,order,ws)
d_target = (I(n^(order+1)) + kron(s',kron(kron(s',s'),t)))\d_orig
@test d ≈ d_target


n1 = 4
n2 = 3
a = randn(n1,n1)
b = randn(n1,n1)
c = randn(n2,n2)
order = 2
ws = EyePlusAtKronBWs(n1,n1,n2,order)
d = randn(n1,n2^order)
a_orig = copy(a)
b_orig = copy(b)
c_orig = copy(c)
d_orig = copy(d)

generalized_sylvester_solver!(a,b,c,vec(d),2,ws)
@test a_orig*d + b_orig*d*kron(c_orig,c_orig) ≈ d_orig
@test d ≈ reshape((kron(I(n2^order),a_orig) + kron(kron(c_orig',c_orig'),b_orig))\vec(d_orig),n1,n2^order)
end

function f(t,s,d,order,ws)
    for i = 1:100
        solver!(t,s,d,order,ws)
    end
end
    
#@profile  f(t,s,d,order,ws)
#Profile.print(combine=true,sortedby=:count)

function kron_power(x,order)
    if order == 0
        return 1
    elseif order == 1
        return x
    else
        m, n = size(x)
        y = Matrix{Float64}(undef, m^order, n^order)
        y1 = similar(y)
        v = view(y, 1:m, 1:n)
        v1 = view(y1, 1:m, 1:n)
        v .= x
        for i = 1:(order-1)
            tmp = v1
            v1 = v
            v = tmp
            v = view(v.parent, 1:m^(i + 1), 1:n^(i + 1))
            kron!(v, v1, x)
        end
    end
    return v.parent
end

function test1(m, n, order)

    a = randn(m,m)
    b = randn(m,m)
    c = randn(n,n)
    d = randn(m*n^order)
    alpha = randn()
    beta = randn()
    t = lu(a\b).U
    s = lu(c).U
    t2 = QuasiUpperTriangular(t*t)
    s2 = QuasiUpperTriangular(s*s)
    t = QuasiUpperTriangular(t)
    s = QuasiUpperTriangular(s)
    r = 1.0
    nd = m*n^(order-1)
    
    ws = GeneralizedSylvesterWs(m, m, n, order)

    @test kron_power(a,1) == a
    @test kron_power(a,3) ≈ kron(kron(a,a),a) 

    println("SOLVI_REAL_ELIMINATE")
    drange = 1:nd
    for index = 1:(n-1)
        d_orig = copy(d)
        @time solvi_real_eliminate!(index, n, nd, drange, order-1, r, t, s, d, ws)
        d_target = d_orig[(index*nd + 1): m*n^order]  - r*kron(s[index,(index+1):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))]
        @test d_target ≈ d[(index*nd + 1): m*n^order]
        drange = drange .+ nd
    end

    println("SOLVI_COMPLEX_ELIMINATE")
    drange = 1:nd
    for index = 1:2:(n-2)
        println(index)
        d_orig = copy(d)
        @time solvi_complex_eliminate!(index, n, nd, drange, order-1, r, t, s, d, ws)
        d_target = (d_orig[(index+1)*nd .+ 1: m*n^order]  - r*kron(s[index,(index+2):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))]
                    - r*kron(s[index + 1,(index+2):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd + m*n^(order-1) .+ (1:m*n^(order-1))])
        @test d_target ≈ d[(index+1)*nd .+ 1: m*n^order]
        drange = drange .+ 2*nd
    end

    println("SOLVIIP_REAL_ELIMINATE")
    drange = 1:nd
    for index = 1:(n-1)
        d_orig = copy(d)
        @time solviip_real_eliminate!(index, n, nd, drange, order-1, alpha, beta, t, t2, s, s2, d, ws)
        d_target = (d_orig[(index*nd + 1): m*n^order]
                    - 2*alpha*kron(s[index,(index+1):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))]
                    - (alpha^2 + beta^2)*kron(s2[index,(index+1):n],kron(kron_power(s2',order-1),t2))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))])
        @test d_target ≈ d[(index*nd + 1): m*n^order]
        drange = drange .+ nd
    end

    println("SOLVIIP_COMPLEX_ELIMINATE")
    drange = 1:nd
    for index = 1:2:(n-2)
        d_orig = copy(d)
        @time solviip_complex_eliminate!(index, n, nd, drange, order-1, alpha, beta, t, t2, s, s2, d, ws)
        d_target = (d_orig[(index+1)*nd + 1: m*n^order]
                    - 2*alpha*kron(s[index,(index+2):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))]
                    - (alpha^2 + beta^2)*kron(s2[index,(index+2):n],kron(kron_power(s2',order-1),t2))*d_orig[(index-1)*nd .+ (1:m*n^(order-1))]
                    - 2*alpha*kron(s[index + 1,(index+2):n],kron(kron_power(s',order-1),t))*d_orig[(index-1)*nd + m*n^(order-1) .+ (1:m*n^(order-1))]
                    - (alpha^2 + beta^2)*kron(s2[index + 1,(index+2):n],kron(kron_power(s2',order-1),t2))*d_orig[(index-1)*nd  + m*n^(order-1) .+ (1:m*n^(order-1))])
        @test d_target ≈ d[(index+1)*nd .+ 1: m*n^order]
        drange = drange .+ 2*nd
    end

    println("TRANSFORMATION1")
    a = randn()
    b1 = -rand()
    b2 = rand()
    d_orig = copy(d[1:2*m*n^(order-1)])
    @time transformation1(a,b1,b2,order-1,t,s,d,ws)
    d_target = d_orig + kron([a -b1; -b2 a],kron(kron_power(s',order-1),t))*d_orig
    @test d[1:2*m*n^(order-1)] ≈ d_target

    if order > 1
        println("TRANSFORMATION2")
        # check handling of order and nd
        a = randn()
        b1 = -rand()
        b2 = rand()
        r1 = randn()
        r2 = randn()
        d1 = randn(2*m*n^(order-1))
        nd1 = m*n^(order-2)
        d_orig = copy(d1[1:2*m*n^(order-2)])
        @time transform2(r1, r2, a, b1, b2, nd1, order-1, t, t2, s, s2, d1, ws)
        d_target = (I(2*m*n^(order-2)) + 2*r1*kron([a b1; b2 a],kron(kron_power(s',(order-2)),t))
                    + (r1*r1 + r2*r2)*kron([a b1; b2 a]*[a b1; b2 a],kron(kron_power(s2',(order-2)),t2)))*d_orig
        @test d1[1:2*m*n^(order-2)] ≈ d_target
    end
    
    println("SOLVE1 real eigenvalues case")
    d_orig = copy(d)
    sreal = QuasiUpperTriangular(triu(s))
    ts = triu(s)
    sreal2 = QuasiUpperTriangular(ts*ts)
    @time solve1!(r,order,t,t2,sreal,sreal2,d,ws)
    d_target = (I(m*n^order) + kron(kron_power(sreal',order),t))\d_orig
    @test d ≈ d_target

    if n == 3
        alpha = randn()
        beta1 = rand()
        beta2 = rand()
        d_orig = copy(d)

        println("SOLVEIIP beta == 0 s, with real eigenvalues")
        beta1 = 0.0 
        @time solviip(alpha, beta1, order, t, t2, sreal, sreal2, d, ws)
        d_target = (I(m*n^order) + 2*alpha*kron(kron_power(sreal',order),t) + (alpha*alpha + beta1*beta1)*kron(kron_power(sreal2',order),t2))\d_orig
        @test d ≈ d_target
        
        scplx = copy(s)
        scplx[1:2,1:2] = [alpha -beta1; beta2 alpha]
        scplx[2,3] = 0
        scplx2 = QuasiUpperTriangular(scplx*scplx)
        
        println("SOLVE1 complex eigenvalues case")
        d = randn(m*n^order)
        d_orig = copy(d)
        @time solve1!(r,order,t,t2,scplx,scplx2,d,ws)
        d_target = (I(m*n^order) + kron(kron_power(scplx',order),t))\d_orig
        @test d ≈ d_target
    end
    
    println("SOLVEIIP with real eigenvalue")
    d_orig = copy(d)
    @time solviip(alpha, beta, order,
                                       t, t2, sreal, sreal2,
                                       d, ws)
d_target = (I(m*n^order)
            + 2*alpha*kron(kron_power(sreal',order),t)
            + (alpha*alpha + beta*beta)*kron(kron_power(sreal2',order),t2))\d_orig
@test d ≈ d_target

println("SOLVEIIP2")
d = randn(2*m*n^(order-1))
d_orig = copy(d)
a = randn()
b1 = -rand()
b2 = rand()
G = [a b1; b2 a]
@time solviip2(alpha,beta,a,b2,b1,order,t,t2,s,s2,d,ws)
d_target = (I(2*m*n^(order-1)) + 2*alpha*kron(kron(G',kron_power(s',order-1)),t)
            + (alpha*alpha + beta*beta)*kron(kron(G'*G',kron_power(s2',order-1)),t2))\d_orig
@test d ≈ d_target


println("SOLVE1 general case")
d = randn(m*n^order)
d_orig = copy(d)
@time solve1!(r,order,t,t2,s,s2,d,ws)
d_target = (I(m*n^order) + kron(kron_power(s',order),t))\d_orig
@test d ≈ d_target

end

n=6
test1(n,n,2)
test1(n,n,2)
