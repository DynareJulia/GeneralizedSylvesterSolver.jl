import GeneralizedSylvesterSolver: GeneralizedSylvesterWs, solvi_real_eliminate!,
    QuasiUpperTriangular, solve1!, transformation1, transform2,
    generalized_sylvester_solver!, solvi_complex_eliminate!, solviip_complex_eliminate!, solver!,
    solvii, solviip, solviip2, solviip_real_eliminate!

using LinearAlgebra    
using Random
using Test

    
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

sreal = [2.0 1.0; 0.0 -3.0]
scplx = [-2.0 3.0; -1.0 -2.0]
@assert scplx[1,1] == scplx[2,2] "not a Schur block"
@assert !(scplx[2,1]*scplx[1,2] > 0) "not a Schur block"

treal = [1.0 3.0; 0.0 -2.0]
tcplx = [1.0 3.0; 1.0 1.0]
o2 = zeros(2,2)

matrices= [ (sreal, treal),
            (scplx, treal),
            (scplx, tcplx),
            (scplx, tcplx),
            ([sreal scplx; o2 sreal], [treal treal; o2 treal]),
            ([sreal scplx; o2 scplx], [treal treal; o2 treal]),
            ([scplx scplx; o2 sreal], [treal treal; o2 treal]),
            ([sreal scplx; o2 sreal], [treal tcplx; o2 treal]),
            ([sreal scplx; o2 scplx], [treal tcplx; o2 treal]),
            ([scplx scplx; o2 sreal], [treal tcplx; o2 treal]),
            ([sreal scplx; o2 sreal], [tcplx tcplx; o2 tcplx]),
            ([sreal scplx; o2 scplx], [tcplx tcplx; o2 tcplx]),
            ([scplx scplx; o2 sreal], [tcplx tcplx; o2 tcplx])
           ]


@testset "Generalized Sylvester Equation $depth, $mat " for depth = 0:4, mat in matrices
    s = QuasiUpperTriangular(mat[1])
    t = QuasiUpperTriangular(mat[2])
    s2 = QuasiUpperTriangular(s*s)
    t2 = QuasiUpperTriangular(t*t)
    n = size(s,1)
    m = size(t,1)
    d_orig = randn(m*n^depth)
    r = 1.5
        
    ws = GeneralizedSylvesterWs(m, m, n, depth)

    @testset "kron_power" begin
        @test kron_power(s,1) == s
        @test kron_power(s,3) ≈ kron(kron(s,s),s) 
    end

    if depth > 0
        nd = m*n^(depth-1)
        @testset "SOLVI_REAL_ELIMINATE" begin
            drange = 1:nd
            d = copy(d_orig)
            for index = 1:(n-1)
                d_target = d[(index*nd + 1): m*n^depth]  - r*kron(s[index,(index+1):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd .+ (1:m*n^(depth-1))]
                solvi_real_eliminate!(index, n, nd, drange, depth-1, r, t, s, d, ws)
                @test d_target ≈ d[(index*nd + 1): m*n^depth]
                drange = drange .+ nd
            end
        end

        @testset "SOLVI_COMPLEX_ELIMINATE" begin
            drange = 1:nd
            d = copy(d_orig)
            for index = 1:2:(n-2)
                d_target = (d[(index+1)*nd .+ 1: m*n^depth]  - r*kron(s[index,(index+2):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd .+ (1:m*n^(depth-1))]
                            - r*kron(s[index + 1,(index+2):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd + m*n^(depth-1) .+ (1:m*n^(depth-1))])
                solvi_complex_eliminate!(index, n, nd, drange, depth-1, r, t, s, d, ws)
                @test d_target ≈ d[(index+1)*nd .+ 1: m*n^depth]
                drange = drange .+ 2*nd
            end
        end

        @testset "SOLVIIP_REAL_ELIMINATE" begin 
            drange = 1:nd
            d = copy(d_orig)
            alpha = randn()
            beta = randn()
            for index = 1:(n-1)
                d_target = (d[(index*nd + 1): m*n^depth]
                            - 2*alpha*kron(s[index,(index+1):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd .+ (1:m*n^(depth-1))]
                            - (alpha^2 + beta^2)*kron(s2[index,(index+1):n],kron(kron_power(s2',depth-1),t2))*d[(index-1)*nd .+ (1:m*n^(depth-1))])
                solviip_real_eliminate!(index, n, nd, drange, depth-1, alpha, beta, t, t2, s, s2, d, ws)
                @test d_target ≈ d[(index*nd + 1): m*n^depth]
                drange = drange .+ nd
            end
        end

        @testset "SOLVIIP_COMPLEX_ELIMINATE" begin
            drange = 1:nd
            d = copy(d_orig)
            alpha = randn()
            beta = randn()
            for index = 1:2:(n-2)
                d_target = (d[(index+1)*nd + 1: m*n^depth]
                            - 2*alpha*kron(s[index,(index+2):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd .+ (1:m*n^(depth-1))]
                            - (alpha^2 + beta^2)*kron(s2[index,(index+2):n],kron(kron_power(s2',depth-1),t2))*d[(index-1)*nd .+ (1:m*n^(depth-1))]
                            - 2*alpha*kron(s[index + 1,(index+2):n],kron(kron_power(s',depth-1),t))*d[(index-1)*nd + m*n^(depth-1) .+ (1:m*n^(depth-1))]
                            - (alpha^2 + beta^2)*kron(s2[index + 1,(index+2):n],kron(kron_power(s2',depth-1),t2))*d[(index-1)*nd  + m*n^(depth-1) .+ (1:m*n^(depth-1))])
                solviip_complex_eliminate!(index, n, nd, drange, depth-1, alpha, beta, t, t2, s, s2, d, ws)
                @test d_target ≈ d[(index+1)*nd .+ 1: m*n^depth]
                drange = drange .+ 2*nd
            end
        end

        @testset "SOLVEIIP2" begin
            alpha = randn()
            beta = randn()
            a = randn()
            b1 = -rand()
            b2 = rand()
            G = [a b1; b2 a]
            nd1 = 2*m*n^(depth-1)
            d = copy(d_orig[1:nd1])
            d_target = (I(nd1) + 2*alpha*kron(kron(G',kron_power(s',depth-1)),t)
                        + (alpha*alpha + beta*beta)*kron(kron(G'*G',kron_power(s2',depth-1)),t2))\d
            solviip2(alpha, beta, a, b2, b1, depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
        end
        
        @testset "TRANSFORMATION1" begin
            a = randn()
            b1 = -rand()
            b2 = rand()
            nd1 = 2*m*n^(depth -1)
            d = copy(d_orig[1:nd1])
            d_target = d + kron([a -b1; -b2 a], kron(kron_power(s', depth - 1), t))*d
            transformation1(a, b1, b2, depth - 1, t, s, d, ws)
            @test d ≈ d_target
        end

        @testset "TRANSFORMATION2" begin
            a = randn()
            b1 = -rand()
            b2 = rand()
            r1 = randn()
            r2 = randn()
            nd1 = 2*m*n^(depth -1)
            d = copy(d_orig[1:nd1])
            d_target = (I(nd1) + 2*r1*kron([a b1; b2 a],kron(kron_power(s', depth -1), t))
                        + (r1*r1 + r2*r2)*kron([a b1; b2 a]*[a b1; b2 a], kron(kron_power(s2', depth - 1), t2)))*d
            transform2(r1, r2, a, b1, b2, m*n^(depth - 1), depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
        end

        @testset "SOLVEII" begin
            alpha = randn()
            beta1 = -rand()
            beta2 = rand()
            G = [alpha beta1; beta2 alpha]
            nd1 = 2*m*n^(depth -1)
            d = copy(d_orig[1:nd1])
            d_target = (I(nd1) + kron(kron(G, kron_power(s',depth-1)), t))\d
            solvii(alpha, beta1, beta2, depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
            beta1 = rand()
            beta2 = -rand()
            G = [alpha beta1; beta2 alpha]
            d = copy(d_orig[1:nd1])
            d_target = (I(nd1) + kron(kron(G, kron_power(s',depth-1)), t))\d
            solvii(alpha, beta1, beta2, depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
        end
        
        @testset "SOLVEIIP beta == 0.0" begin
            alpha = randn()
            beta1 = 0.0
            nd = m*n^(depth - 1)
            d = copy(d_orig[1:nd])
            d_target = (I(nd) + 2*alpha*kron(kron_power(s',depth - 1),t) + (alpha*alpha + beta1*beta1)*kron(kron_power(s2', depth - 1),t2))\d
            solviip(alpha, beta1, depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
        end        

        @testset "SOLVEIIP beta == 2.0" begin
            alpha = randn()
            beta1 = randn()
            nd = m*n^(depth - 1)
            d = copy(d_orig[1:nd])
            d_target = (I(nd) + 2*alpha*kron(kron_power(s',depth - 1),t) + (alpha*alpha + beta1*beta1)*kron(kron_power(s2', depth - 1),t2))\d
            solviip(alpha, beta1, depth - 1, t, t2, s, s2, d, ws)
            @test d ≈ d_target
        end        
    end
    
    @testset "SOLVE1" begin
        d = copy(d_orig)
        r = 1.0
        d_target = (I(m*n^depth) + r*kron(kron_power(s', depth), t))\d
        solve1!(r, depth, t, t2, s, s2, d, ws)
        @test d ≈ d_target
    end

end


n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b).T)
s = QuasiUpperTriangular(schur(c).T)
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
depth = 1
d_orig = randn(n^(depth+1))
d = copy(d_orig)
ws = GeneralizedSylvesterWs(n,n,n,depth)
solver!(t,s,d,depth,ws)
d_target = (I(n^(depth+1)) + kron(s',t))\d_orig
@test d ≈ d_target

n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b).T)
s = QuasiUpperTriangular(schur(c).T)
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
depth = 3
d_orig = randn(n^(depth+1))
d = copy(d_orig)
ws = GeneralizedSylvesterWs(n,n,n,depth)
solver!(t,s,d,depth,ws)
d_target = (I(n^(depth+1)) + kron(s',kron(kron(s',s'),t)))\d_orig
@test d ≈ d_target

n = 4
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)
t = QuasiUpperTriangular(schur(a\b).T)
s = QuasiUpperTriangular(schur(c).T)
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)
depth = 3
ws = GeneralizedSylvesterWs(n,n,n,depth)
d_orig = randn(n^(depth+1))
d = copy(d_orig)
solver!(t,s,d,depth,ws)
d_target = (I(n^(depth+1)) + kron(s',kron(kron(s',s'),t)))\d_orig
@test d ≈ d_target


n1 = 4
n2 = 3
a_orig = randn(n1,n1)
b_orig = randn(n1,n1)
c_orig = randn(n2,n2)

depth = 1
ws = GeneralizedSylvesterWs(n1,n1,n2,depth)
d_orig = randn(n1,n2^depth)
a = copy(a_orig)
b = copy(b_orig)
c = copy(c_orig)
d = copy(d_orig)

d = reshape(d, 4, 3)
generalized_sylvester_solver!(a,b,c,d,1,ws)
@test a_orig*d + b_orig*d*c_orig ≈ d_orig
@test d ≈ reshape((kron(I(n2^depth),a_orig) + kron(c_orig',b_orig))\vec(d_orig),n1,n2^depth)

depth = 2
ws = GeneralizedSylvesterWs(n1,n1,n2,depth)
d_orig = randn(n1,n2^depth)
a = copy(a_orig)
b = copy(b_orig)
c = copy(c_orig)
d = copy(d_orig)

d = reshape(d, 4, 9)
generalized_sylvester_solver!(a,b,c,d,2,ws)
@test a_orig*d + b_orig*d*kron(c_orig,c_orig) ≈ d_orig
@test d ≈ reshape((kron(I(n2^depth),a_orig) + kron(kron(c_orig',c_orig'),b_orig))\vec(d_orig),n1,n2^depth)

function f(t,s,d,depth,ws)
    for i = 1:100
        solver!(t,s,d,depth,ws)
    end
end
    
#@profile  f(t,s,d,depth,ws)
#Profile.print(combine=true,sortedby=:count)

