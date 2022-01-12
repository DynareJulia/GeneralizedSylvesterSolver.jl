import SolveIPlusOrMinusAkronB: IPlusAtKronBWS, solvi_real_eliminate!,
                                QuasiUpperTriangular, solve1!, transformation1,
                                generalized_sylvester_solver!
using Test

n = 3
x = randn(n,n*n)
a = randn(n,n)
b = randn(n,n)
c = randn(n,n)

ws = IPlusAtKronBWS(n, n, n, 2)

index = [n,n]
j = 1
depth = 2
t = QuasiUpperTriangular(lu(a\b).U)
s = QuasiUpperTriangular(lu(c).U)
d = zeros(n^(depth+1))
d[1:n*n] = ones(n*n)
d_orig = copy(d)
nd = n*n^(depth-1)
drange1 = 1:nd
r = 1.0
solvi_real_eliminate!(1, n, nd, drange1, depth - 1, r, t, s, d, ws)
@test ws.work1[drange1] ≈ kron(s', t)*d_orig[drange1]
kst = kron(s', t)
d_target = -kron(s[1, 2:n], kst*d_orig[drange1])
@test d[nd + 1:end] ≈ d_target

index = [1,2]
drange1 = drange1 .+ nd
d = zeros(n^(depth+1))
d[drange1] = ones(nd)
d_orig = copy(d)
r = 1.0

@show drange1
solvi_real_eliminate!(2, n, nd, drange1, depth - 1, r, t, s, d, ws)
d_target = -kron(s[2, 3:n], kst*d_orig[drange1])
@show ws.work1
@show d
@test d[2*nd + 1:end] ≈ d_target


aa = a\b
t = lu(aa).U
s = lu(c).U

r = 1.0
d = randn(n)
d_orig = copy(d)
depth = 0
s2 = QuasiUpperTriangular(s*s)
t2 = QuasiUpperTriangular(t*t)

t = QuasiUpperTriangular(t)
s = QuasiUpperTriangular(s)

solve1!(r, depth,t, t2, s, s2, d, ws)
d_target = (I(n) + t)\d_orig
@test d_target ≈ d

r = 1.0
d = randn(n*n)
d_orig = copy(d)
depth = 1
index = zeros(Int64,depth)
d_block_number = n^depth
solve1!(r, depth,t, t2, s, s2, d, ws)

d_target = (I(n^2) + kron(s,t))\d_orig

@test d ≈ d_target

r = 1.0
d = randn(n*n*n)
d_orig = copy(d)
depth = 2
index = zeros(Int64,depth)
d_block_number = n^depth
solve1!(r, depth,t, t2, s, s2, d, ws)

d_target = (I(n^3) + kron(kron(s,s),t))\d_orig

@test d ≈ d_target

a = 0.5
b1 = 0.1
b2 = 1.3
depth = 2
tt = QuasiUpperTriangular(t)
ss = QuasiUpperTriangular(s)
d = randn(2*n^depth)
d_orig = copy(d)
transformation1(a, b1, b2, depth, tt, ss, d, ws)
d_target = d_orig + kron([a -b1; -b2 a],kron(s,t))*d_orig

@test d ≈ d_target


d = a*x + b*x*kron(c,c)
generalized_sylvester_solver!(a, b, c, d, 2, ws)

@test d ≈ x
