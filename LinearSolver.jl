function LinSolver(k, l, U)
    e = [[1., 0.] [0., 1.]]
    function unitStrain(k,l)
        return e[k,:]*e[l,:]'
    end

    Chia, Chib = unitStrain(k,l)*L*e[1,:], unitStrain(k,l)*L*e[2,:]

    Chicorner = vcat(zeros(2, 1), Chia, Chia+Chib, Chib)
    Chiinternal = zeros(Float64, 2*Ns[2])
    Chimaster = zeros(Float64, 2*Ns[3])
    Chislave = Chimaster + D*Chicorner
    Chi = vcat(Chicorner, Chiinternal, Chimaster, Chislave)

    function Assembly(U::Matrix{Float64})
        gps = gaussPoints2D(2)
        K, Fint = zeros(Float64, Ns[6], Ns[6]), zeros(Float64, Ns[6], 1) 
        for i = 1:Ns[5]
            Ce = vec(DOFs[i, :, :]')
            elK, elFint = zeros(Float64, 8, 8), zeros(Float64, 8, 1)
            for gp in eachrow(gps)
                g, h = gp[2:end]
                B, detJ, PK1_1D, A_2D = constitutiveEqs(U[Ce, :], elCoords[i, :, :], g, h, E, nu)
                elK += gp[1]*(B'*A_2D*B)*detJ
            end
            K[Ce, Ce] += elK
        end
        return sparse(K)
    end

    function solve(K::SparseMatrixCSC, Chicorner)
        Kic = K[2*Ns[1]+1:2*Ns[7], 1:2*Ns[1]]
        Kii = K[2*Ns[1]+1:2*Ns[7], 2*Ns[1]+1:2*Ns[7]]
        Kim = K[2*Ns[1]+1:2*Ns[7], 2*Ns[7]+1:2*Ns[9]]
        Kis = K[2*Ns[1]+1:2*Ns[7], 2*Ns[9]+1:2*Ns[11]]

        Kmc = K[2*Ns[7]+1:2*Ns[9], 1:2*Ns[1]]
        Kmi = K[2*Ns[7]+1:2*Ns[9], 2*Ns[1]+1:2*Ns[7]]
        Kmm = K[2*Ns[7]+1:2*Ns[9], 2*Ns[7]+1:2*Ns[9]]
        Kms = K[2*Ns[7]+1:2*Ns[9], 2*Ns[9]+1:2*Ns[11]]

        Ksc = K[2*Ns[9]+1:2*Ns[11], 1:2*Ns[1]]
        Ksi = K[2*Ns[9]+1:2*Ns[11], 2*Ns[1]+1:2*Ns[7]]
        Ksm = K[2*Ns[9]+1:2*Ns[11], 2*Ns[7]+1:2*Ns[9]]
        Kss = K[2*Ns[9]+1:2*Ns[11], 2*Ns[9]+1:2*Ns[11]]

        Qi = -(Kic+Kis*D)*Chicorner
        Qm = -(Kmc+Kms*D+Ksc+Kss*D)*Chicorner
        Kupper = hcat(Kii, Kim + Kis)
        Klower = hcat(Kmi + Ksi, Kmm + Kms + Ksm + Kss)

        A = vcat(Kupper, Klower)
        Q = vcat(Qi, Qm)

        sol = A \ Q

        return sol[1:2*Ns[2],:], sol[2*Ns[2]+1:end,:]
    end

    K = Assembly(U)
    Chiinternal, Chimaster = solve(K, Chicorner)
    Chislave = Chimaster + D*Chicorner
    Chi = vcat(Chicorner, Chiinternal, Chimaster, Chislave)
    
    return Chi
end