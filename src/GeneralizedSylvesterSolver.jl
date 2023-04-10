module GeneralizedSylvesterSolver
###
# solving a x + b x (c ⊗ c ⊗ ... ⊗ c) = d
# using (I + s^T ⊗ s^T ⊗ ... \otimes s^T \otimes t)x = d
###

using QuasiTriangular
using FastLapackInterface
using KroneckerTools
using LinearAlgebra

export GeneralizedSylvesterWs, generalized_sylvester_solver!

struct GeneralizedSylvesterWs
    ma::Int64
    mb::Int64
    b1::Matrix{Float64}
    c1::Matrix{Float64}
    vs_b::Matrix{Float64}
    vs_c::Matrix{Float64}
    s2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    t2::QuasiUpperTriangular{Float64,Matrix{Float64}}
    work1::Vector{Float64}
    work2::Vector{Float64}
    work3::Vector{Float64}
    result::Matrix{Float64}
<<<<<<< Updated upstream
    linsolve::LUWs
    schur_b::SchurWs
    schur_c::SchurWs
=======
    linsolve_ws::LUWs
    dgees_ws_b::SchurWs
    dgees_ws_c::SchurWs
>>>>>>> Stashed changes
    function GeneralizedSylvesterWs(ma::Int64, mb::Int64, mc::Int64, order::Int64)
        if mb != ma
            DimensionMismatch("a has $ma rows but b has $mb rows")
        end
        b1 = Matrix{Float64}(undef, mb,mb)
        c1 = Matrix{Float64}(undef, mc,mc)
        vs_b = Matrix{Float64}(undef, mb,mb)
        vs_c = Matrix{Float64}(undef, mc,mc)
        s2 = QuasiUpperTriangular(Matrix{Float64}(undef, mc,mc))
        t2 = QuasiUpperTriangular(Matrix{Float64}(undef, mb,mb))
<<<<<<< Updated upstream
        linsolve = LUWs(ma)
        schur_b = SchurWs(b1)
        schur_c = SchurWs(c1)
=======
        linsolve_ws = LUWs(ma)
        dgees_ws_b = SchurWs(vs_b)
        dgees_ws_c = SchurWs(vs_c)
>>>>>>> Stashed changes
        work1 = Vector{Float64}(undef, ma*mc^order)
        work2 = Vector{Float64}(undef, ma*mc^order)
        work3 = Vector{Float64}(undef, ma*mc^order)
        result = Matrix{Float64}(undef, ma,mc^order)
        new(ma, mb, b1, c1, vs_b, vs_c, s2, t2, work1, work2, work3, result, linsolve, schur_b, schur_c)
    end
end

function generalized_sylvester_solver!(a::AbstractMatrix,b::AbstractMatrix,c::AbstractMatrix,
                                   d::AbstractMatrix,order::Int64,ws::GeneralizedSylvesterWs)
    copy!(ws.b1,b)
    copy!(ws.c1,c)
<<<<<<< Updated upstream
    ws_a = copy(a) # avoid mutating inputs

    # linsolve_core!(a, ws.b1, ws.linsolve_ws)
    factors = LinearAlgebra.LU( LAPACK.getrf!(ws.linsolve, ws_a)... )
    ldiv!(factors, ws.b1) #confirmed

    # linsolve_core_no_lu!(a, d, ws.linsolve_ws)
    factors = LinearAlgebra.LU( LAPACK.getrf!(ws.linsolve, ws_a)... )
    ldiv!(factors, d) #confirmed
    
    # dgees!(ws.dgees_ws_b, ws.b1)
    Schur(LAPACK.gees!(ws.schur_b, 'V', ws.b1)...) #confirmed
    # dgees!(ws.dgees_ws_c, ws.c1)
    Schur(LAPACK.gees!(ws.schur_c, 'V', ws.c1)...) #confirmed

=======
    linsolve_core!(a, ws.b1, ws.linsolve_ws)
    linsolve_core_no_lu!(a, d, ws.linsolve_ws)
    gees!(ws.dgees_ws_b,ws.b1)
    gees!(ws.dgees_ws_c,ws.c1)
>>>>>>> Stashed changes
    t = QuasiUpperTriangular(ws.b1)
    mul!(ws.t2,t,t) #confirmed
    s = QuasiUpperTriangular(ws.c1)
<<<<<<< Updated upstream
    mul!(ws.s2,s,s) #confirmed
    #confirmed as  ws.result = ws.dgees_ws_b.vs' * d * kron(ws.dgees_ws_c.vs, ws.dgees_ws_c.vs)
    at_mul_b_kron_c!(ws.result, ws.schur_b.vs, d, ws.schur_c.vs, order, ws.work2, ws.work3)
    copy!(d, ws.result)

=======
    A_mul_B!(ws.s2,s,s)
    at_mul_b_kron_c!(ws.result, ws.dgees_ws_b.vs, d, ws.dgees_ws_c.vs, order, ws.work2, ws.work3)
    copy!(d,ws.result)
    @show d
>>>>>>> Stashed changes
    solve1!(1.0, order, t, ws.t2, s, ws.s2, vec(d), ws)
    a_mul_b_kron_ct!(ws.result, ws.schur_b.vs, d, ws.schur_c.vs, order, ws.work2, ws.work3)
    copy!(d, reshape(ws.result, size(a, 1), size(c, 2)^order))
end


function solver!(t::QuasiUpperTriangular,s::QuasiUpperTriangular,d::AbstractVector,order::Int64,ws::GeneralizedSylvesterWs)
    s2 = QuasiUpperTriangular(s*s)
    t2 = QuasiUpperTriangular(t*t)
    solve1!(1.0,order,t,t2,s,s2,d,ws)
    d
end

function solve1!(r::Float64, depth::Int64, t::AbstractArray{Float64,2}, t2::AbstractArray{Float64,2}, s::AbstractArray{Float64,2}, s2::AbstractArray{Float64,2}, d::AbstractVector{Float64}, ws::GeneralizedSylvesterWs)
    m = size(t,2)
    n = size(s,1)
    if depth == 0
        I_plus_rA_ldiv_B!(r,t,d)
    else
        nd = m*n^(depth-1)
        nd2 = 2*nd
        drange1 = 1:nd
        drange2 = 1:nd2
        i = 1
        while i <= n
            if i == n || s[i+1,i] == 0
                dv = view(d,drange1)
                solve1!(r*s[i,i],depth-1,t,t2,s,s2,dv,ws)
                if i < n
                    solvi_real_eliminate!(i,n,nd,drange1,depth-1,r,t,s,d,ws)
                end
                drange1 = drange1 .+ nd
                drange2 = drange2 .+ nd
                i += 1
            else
                dv = view(d,drange2)
                solvii(r*s[i,i],r*s[i+1,i],r*s[i,i+1],depth-1,t,t2,s,s2,dv,ws)
                if i < n - 1
                    solvi_complex_eliminate!(i,n,nd,drange1,depth-1,r,t,s,d,ws)
                end
                drange1 = drange1 .+ nd2
                drange2 = drange2 .+ nd2
                i += 2
            end
        end
    end
end

"""
solvi_real_eliminate!(i::Int64, n::Int64, nd::Int64, drange::UnitRange{Int64},
                      depth::Int64, r::Float64, t::QuasiUpperTriangular, s::QuasiUpperTriangular,
                      d::AbstractVector, ws::GeneralizedSylvesterWs)
    updates d[k] with d[k] - ∑_{j=i+1}^n r*s[i, j]*(s ⊗ s ⊗ ... ⊗ s ⊗ t)*d with
    depth occurences of s
"""
function solvi_real_eliminate!(i::Int64, n::Int64, nd::Int64, drange::UnitRange{Int64},
                               depth::Int64, r::Float64, t::QuasiUpperTriangular, s::QuasiUpperTriangular,
                               d::AbstractVector, ws::GeneralizedSylvesterWs)
    work1 = ws.work1
    work2 = ws.work2
    work3 = ws.work3
    kron_at_kron_b_mul_c!(work1,1,s,depth,t,d,drange[1],work2,work3,1)
    k1 = drange[1] + nd
    @inbounds for j = i+1:n
        m = r*s[i,j]
        @simd for k2 = 1:nd
            d[k1] -= m*work1[k2]
            k1 += 1
        end
    end
end

function solvi_complex_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                                  depth::Int64,r::Float64,t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                                  d::AbstractVector,ws::GeneralizedSylvesterWs)
    work1 = ws.work1
    work2 = ws.work2
    work3 = ws.work3
    kron_at_kron_b_mul_c!(work1, 1, s, depth, t, d, drange[1], work2, work3, 1)
    drange = drange .+ nd
    kron_at_kron_b_mul_c!(work1, nd+1, s, depth, t, d, drange[1], work2, work3, 1)
    k1 = drange[1] + nd
    @inbounds for j = i + 2 : n
        m1 = r*s[i,j]
        m2 = r*s[i+1,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*work1[k2] + m2*work1[k2 + nd]
            k1 += 1
        end
    end 
end

function solvii(alpha::Float64,beta1::Float64,beta2::Float64,depth::Int64,
                t::QuasiUpperTriangular,t2::QuasiUpperTriangular,s::QuasiUpperTriangular,
                s2::QuasiUpperTriangular,d::AbstractVector,ws::GeneralizedSylvesterWs)
    m = size(t,2)
    n = size(s,1)
    nd = m*n^depth
    transformation1(alpha,beta1,beta2,depth,t,s,d,ws)
    dv = view(d,1:nd)
    solviip(alpha,sqrt(-beta1*beta2),depth,t,t2,s,s2,dv,ws)
    dv = view(d,nd+1:2*nd)
    solviip(alpha,sqrt(-beta1*beta2),depth,t,t2,s,s2,dv,ws)
end

function transformation1(a::Float64,b1::Float64,b2::Float64,depth::Int64,
                         t::QuasiUpperTriangular,s::QuasiUpperTriangular,
                         d::AbstractVector,ws::GeneralizedSylvesterWs)
    m = size(t, 2)
    n = size(s, 1)
    nd = m*n^depth
    copyto!(ws.work3, d)
    drange = 1:nd
    d1 = view(ws.work3, drange)
    d2 = view(ws.work3, drange .+ nd)
    work = view(ws.work2,drange)
    kron_at_kron_b_mul_c!(s,depth,t,d1,work)
    kron_at_kron_b_mul_c!(s,depth,t,d2,work)
    @inbounds @simd for i = drange
        d[i] += a*d1[i] - b1*d2[i]
        d[i+nd] += -b2*d1[i] + a*d2[i] 
    end
end

diag_zero_sq = 1e-30

function solviip(alpha::Float64,beta::Float64,depth::Int64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                 s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::GeneralizedSylvesterWs)
    m = size(t,2)
    n = size(s,1)
    if beta*beta < diag_zero_sq
        @show "OK100"
        solve1!(alpha,depth,t,t2,s,s2,d,ws)
        solve1!(alpha,depth,t,t2,s,s2,d,ws)
        return
    end

    if depth == 0
        @show "OK0"
        @show size(t), size(d)
        d_orig = copy(d)
        @show d
        I_plus_rA_plus_sB_ldiv_C!(2*alpha,alpha*alpha+beta*beta,t,t2,d)
        @show d
        @show d_orig
        nt = size(t, 1)
        @show d[1:nt] - (I(nt) + 2*alpha*t + (alpha*alpha+beta*beta)*t2)\d_orig[1:n]
        @show d
    else
        nd = m*n^(depth-1)
        nd2 = 2*nd
        drange1 = 1:nd
        drange2 = 1:nd2
        i = 1
        while i <= n
            if i == n || s[i+1,i] == 0
                @show "OK1"
                dv = view(d,drange1)
                if s[i,i]*s[i,i]*(alpha*alpha+beta*beta) > diag_zero_sq
                    solviip(s[i,i]*alpha,s[i,i]*beta,depth-1,t,t2,s,s2,dv,ws)
                end
                if i < n
                    solviip_real_eliminate!(i,n,nd,drange1,depth-1,alpha,beta,t,t2,s,s2,d,ws)
                end
                drange1 = drange1 .+ nd
                drange2 = drange2 .+ nd
                i += 1
            else
                @show "OK2"
                dv = view(d,drange2)
                solviip2(alpha,beta,s[i,i],s[i+1,i],s[i,i+1],depth,t,t2,s,s2,dv,ws)
                if i < n - 1
                    solviip_complex_eliminate!(i,n,nd,drange1,depth-1,alpha,beta,t,t2,s,s2,d,ws)
                end
                drange1 = drange1 .+ nd2
                drange2 = drange2 .+ nd2
                i += 2
            end
        end
    end
end

function solviip_real_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},
                                 depth::Int64,alpha::Float64,beta::Float64,t::QuasiUpperTriangular,
                                 t2::QuasiUpperTriangular,s::QuasiUpperTriangular,
                                 s2::QuasiUpperTriangular,d::AbstractVector,ws::GeneralizedSylvesterWs)
    y1 = view(ws.work1,drange)
    y2 = view(ws.work1,drange .+ nd)
    copyto!(y1,1,d,drange[1],nd)
    copyto!(y2,1,d,drange[1],nd)
    work = view(ws.work2,1:length(drange))
    kron_at_kron_b_mul_c!(s,depth,t,y1,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y2,work)
    k1 = drange[1] + nd
    @inbounds for j = i+1:n
        m1 = 2*alpha*s[i,j]
        m2 = (alpha*alpha+beta*beta)*s2[i,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*y1[k2] + m2*y2[k2]
            k1 += 1
        end
    end 
end

function solviip2(alpha::Float64,beta::Float64,gamma::Float64,delta1::Float64,delta2::Float64,
                  depth::Int64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                  s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::GeneralizedSylvesterWs)
    m = size(t,2)
    n = size(s,1)
    aspds = alpha*alpha + beta*beta
    gspds = gamma*gamma - delta1*delta2
    nd = m*n^(depth-1)
    dv1 = view(d,1:nd)
    dv2 = view(d,nd+1:2*nd)
    if aspds*gspds > diag_zero_sq
        transform2(alpha, beta, gamma, -delta1, -delta2, nd, depth, t, t2, s, s2, d,ws)

        delta = sqrt(-delta1*delta2)
	a1 = alpha*gamma - beta*delta
	b1 = alpha*delta + gamma*beta
	a2 = alpha*gamma + beta*delta
	b2 = alpha*delta - gamma*beta
	solviip(a2, b2, depth-1, t, t2, s, s2, dv1, ws);
	solviip(a1, b1, depth-1, t, t2, s, s2, dv1, ws);
        solviip(a2, b2, depth-1, t, t2, s, s2, dv2, ws);
        solviip(a1, b1, depth-1, t, t2, s, s2, dv2, ws);
    end
end

function transform2(alpha::Float64, beta::Float64, gamma::Float64, delta1::Float64, delta2::Float64,
                    nd::Int64, depth::Int64, t::QuasiUpperTriangular, t2::QuasiUpperTriangular,
                    s::QuasiUpperTriangular, s2::QuasiUpperTriangular, d::AbstractVector, ws::GeneralizedSylvesterWs)
    d1 = ws.work1
    kron_at_kron_b_mul_c!(d1,1,s,depth-1,t,d,1,ws.work2,ws.work3,1)
    kron_at_kron_b_mul_c!(d1,nd+1,s,depth-1,t,d,nd+1,ws.work2,ws.work3,1)

    m1 = 2*alpha*gamma
    m2 = 2*alpha*delta1
    @inbounds @simd for i = 1:nd
        dtmp = d1[i]
        d1[i] = d[i] + m1*dtmp + m2*d1[i+nd]
        d1[i+nd] = d[i+nd] + 2*alpha*(delta2*dtmp + gamma*d1[i+nd])
    end

    kron_at_kron_b_mul_c!(d,1,s2,depth-1,t2,d,1,ws.work2,ws.work3,1)
    kron_at_kron_b_mul_c!(d,nd+1,s2,depth-1,t2,d,nd+1,ws.work2,ws.work3,1)
    
    aspds = alpha*alpha + beta*beta;
    gspds = gamma*gamma + delta1*delta2;
    m1 = aspds*gspds
    m2 = 2*aspds*gamma*delta1
    m3 = 2*aspds*gamma*delta2
    @inbounds @simd for i = 1:nd
        dtmp = d[i]
        d[i] = d1[i] + m1*dtmp + m2*d[i+nd]
        d[i+nd] = d1[i+nd] + m3*dtmp + m1*d[i+nd]
    end
end

"""
    solviip_complex_eliminate!(i,n,nd,drange,depth,alpha,beta,t,t2,s,s2,d)

perfoms elimination after solving for a complex diagonal block of size 2*n^depth

d n^(depth+2) x 1

The solution is stored in d[drange; drange + nd]

The function updates d[i*nd+1:n*nd]
"""
function solviip_complex_eliminate!(i::Int64,n::Int64,nd::Int64,drange::UnitRange{Int64},depth::Int64,
                                    alpha::Float64,beta::Float64,t::QuasiUpperTriangular,t2::QuasiUpperTriangular,
                                    s::QuasiUpperTriangular,s2::QuasiUpperTriangular,d::AbstractVector,ws::GeneralizedSylvesterWs)
    y11 = view(ws.work1,drange)
    y12 = view(ws.work1, drange .+ nd)
    copyto!(y11,1,d,drange[1],nd)
    copyto!(y12,1,d,drange[1],nd)
    drange = drange .+ nd
    y21 = view(ws.work2,drange)
    y22 = view(ws.work2, drange .+ nd)
    copyto!(y21,1,d,drange[1],nd)
    copyto!(y22,1,d,drange[1],nd)
    work = view(ws.work3,drange)

    kron_at_kron_b_mul_c!(s,depth,t,y11,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y12,work)
    kron_at_kron_b_mul_c!(s,depth,t,y21,work)
    kron_at_kron_b_mul_c!(s2,depth,t2,y22,work)

    alpha2beta2 = alpha*alpha + beta*beta
    k1 = drange[1] + nd
    @inbounds for j = i+2:n
        m1 = 2*alpha*s[i,j]
        m2 = alpha2beta2*s2[i,j]
        m3 = 2*alpha*s[i+1,j]
        m4 = alpha2beta2*s2[i+1,j]
        @simd for k2 = 1:nd
            d[k1] -= m1*y11[k2] + m2*y12[k2] + m3*y21[k2] + m4*y22[k2]
            k1 += 1
        end
    end 
end

end
