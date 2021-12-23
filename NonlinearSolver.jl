using SparseArrays, LinearAlgebra, Luxor, Colors, ProgressMeter

function NonLinSolver(Fm::Matrix{Float64}, Ns, L::Float64, E::Float64, nu::Float64, DT::Float64)
    L1, L2 = [L, 0.], [0., L]
    U, Fext, Fint = zeros(Float64, Ns[6], 1), zeros(Float64, Ns[6], 1), zeros(Float64, Ns[6], 1)
    gps = gaussPoints2D(2)
    function Assembly(U::Matrix{Float64})
        # Update the elements and assemble the system matrices.
        Fint, Krow, Kcol, K = zeros(Float64, Ns[6], 1), Int64[], Int64[], zeros(Float64, Ns[6], Ns[6])
        for i = 1:Ns[5]
            Ce = vec(DOFs[i, :, :]')
            elFint, elK = Quad4(U[Ce, :], i, elCoords[i, :, :], E, nu)
            Fint[Ce] += elFint
            K[Ce, Ce] += elK
        end
        K = sparse(K)
        Kims = K[2*Ns[1]+1:end, 2*Ns[1]+1:end]
        Rims = (Fint-Fext)[2*Ns[1]+1:end, :]
        return Fint, Kims, Rims
    end
    function generateEq(Kims::SparseMatrixCSC, Rims::Matrix{Float64})
        Kii = Kims[1:2*Ns[2], 1:2*Ns[2]]
        Kim = Kims[1:2*Ns[2], 2*Ns[2]+1:2*Ns[8]]
        Kis = Kims[1:2*Ns[2], 2*Ns[8]+1:2*Ns[10]]

        Kmi = Kims[2*Ns[2]+1:2*Ns[8], 1:2*Ns[2]]
        Kmm = Kims[2*Ns[2]+1:2*Ns[8], 2*Ns[2]+1:2*Ns[8]]
        Kms = Kims[2*Ns[2]+1:2*Ns[8], 2*Ns[8]+1:2*Ns[10]]

        Ksi = Kims[2*Ns[8]+1:2*Ns[10], 1:2*Ns[2]]
        Ksm = Kims[2*Ns[8]+1:2*Ns[10], 2*Ns[2]+1:2*Ns[8]]
        Kss = Kims[2*Ns[8]+1:2*Ns[10], 2*Ns[8]+1:2*Ns[10]]

        Ri = Rims[1:2*Ns[2], :]
        Rm = Rims[2*Ns[2]+1:2*Ns[8], :]
        Rs = Rims[2*Ns[8]+1:2*Ns[10], :]

        Kupper = hcat(Kii, Kim + Kis)
        Klower = hcat(Kmi + Ksi, Kmm + Ksm + Kms + Kss)

        K = vcat(Kupper, Klower)
        R = vcat(Ri, Rm + Rs)
        RNorm = norm(R)

        return K, R, RNorm
    end

    function Iterate(TOL::Float64, DT::Float64, Draw::Bool)
        t, tFinal = (0.0, 1.0)
        A(t) = (t/tFinal)*(Fm-I)[1:2, 1:2]

        p1 = Progress(Int64(tFinal/DT))
        while round(t, digits = 2) < tFinal
            t += DT
            Ua, Ub = A(t)*L1, A(t)*L2
            Ucorner = vcat(zeros(2, 1), Ua, Ua + Ub, Ub)
            Uinternal = U[2*Ns[1]+1:2*Ns[7], :]
            Umaster = U[2*Ns[7]+1:2*Ns[9], :]
            Uslave = Umaster + D * Ucorner
            U = vcat(Ucorner, Uinternal, Umaster, Uslave)
            Fint, Kims, Rims = Assembly(U)
            RNorm = Inf
            #@printf("-----t: %.2f-----\n",t)
            #p2 = ProgressThresh(TOL, "Minimizing: ")
            #update!(p2, RNorm)
            while RNorm > TOL
                K, R, RNorm = generateEq(Kims, Rims)
                dU = K \ -R
                Uinternal += dU[1:2*Ns[2], :]    # dUInternal
                Umaster += dU[2*Ns[2]+1:end, :]  # dUmaster
                Uslave = Umaster + D * Ucorner
                U = vcat(Ucorner, Uinternal, Umaster, Uslave)
                Fint, Kims, Rims = Assembly(U)
                #@printf("  | RNorm: %.3f \n", RNorm)
                #update!(p2, RNorm)
            end
            next!(p1)
        end
        if Draw
            deformed = elCoords + U[DOFs]
            maxUx, maxUy = max(maximum(elCoords[:,:,1]), maximum(deformed[:,:,1])), max(maximum(elCoords[:,:,2]), maximum(deformed[:,:,2]))
            scl = 200
            Drawing(Int(round(maxUx*scl*1.4)), Int(round(maxUy*scl*1.4)))
            background("white")
            origin(maxUx*scl*0.2, maxUy*scl*1.2)
            for i = 1:Ns[5]
                undefPoints = [Point(elCoords[i, j, 1]*scl, -elCoords[i, j, 2]*scl) for j in 1:4]
                defPoints = [Point(deformed[i, j, 1]*scl, -deformed[i, j, 2]*scl) for j in 1:4]
                sethue("black")
                poly(undefPoints, :stroke, close = true)
                sethue("blue")
                poly(defPoints, :stroke, close = true)
            end
            finish();
        end
        return U
    end
    return Iterate(1e-6, DT, true)
end