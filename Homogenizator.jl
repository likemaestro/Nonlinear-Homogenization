function Homogenization(U, Chis, N, E, nu)
    function homogenizeP()
        function MacroP(dir::Int64, U::Matrix{Float64}, elID::Int64, coords::Matrix{Float64}, E::Float64, nu::Float64)
            gps = gaussPoints1D(2)
            if dir == 1
                P = zeros(Float64, 4, 1)
                for gp in eachrow(gps)
                    B, detJ, PK1_1D, A_2D = constitutiveEqs(U, coords, 1/2, gp[2], E, nu)
                    P += gp[1]*PK1_1D*(1/2)
                end
                P11, P21 = P[1], P[4]
                return P11, P21
            elseif dir == 2
                P = zeros(Float64, 4, 1)
                for gp in eachrow(gps)
                    B, detJ, PK1_1D, A_2D = constitutiveEqs(U, coords, gp[2], 1/2, E, nu)
                    P += gp[1]*PK1_1D*(1/2)
                end
                P22, P12 = P[2], P[3]
                return P22, P12
            end
        end
        
        P11, P21 = 0, 0
        for i in N:N:N*N
            Ce = vec(DOFs[i, :, :]')
            P11_i, P21_i = MacroP(1, U[Ce, :], i, elCoords[i, :, :], E, nu)
            P11 += P11_i
            P21 += P21_i
        end
    
        P22, P12 = 0, 0
        for i in N*(N-1)+1:1:N*N
            Ce = vec(DOFs[i, :, :]')
            P22_i, P12_i = MacroP(2, U[Ce, :], i, elCoords[i, :, :], E, nu)
            P22 += P22_i
            P12 += P12_i
        end
        return L/(L*L)*hcat(P11, P22, P12, P21)
    end

    function homogenizeA()
        function MacroA(dir::Int64, U::Matrix{Float64}, elID::Int64, coords::Matrix{Float64}, E::Float64, nu::Float64)
            gps = gaussPoints1D(2)
            if dir == 1
                Ael = zeros(Float64, 4, 4)
                Ce = vec(DOFs[elID, :, :]')
                for i in 1:4
                    for j in 1:4
                        for gp in eachrow(gps)
                            B, detJ, PK1_1D, A_2D = constitutiveEqs(U, coords, 1/2, gp[2], E, nu)
                            Ql = B*Chis[i][Ce,:]
                            Qr = B*Chis[j][Ce,:]
                            Ael[i,j] += gp[1]*(Ql'*A_2D*Qr)[1]*(dx/2)
                        end
                    end
                end
                return Ael
            elseif dir == 2
                Ael = zeros(Float64, 4, 4)
                Ce = vec(DOFs[elID, :, :]')
                for i in 1:4
                    for j in 1:4
                        for gp in eachrow(gps)
                            B, detJ, PK1_1D, A_2D = constitutiveEqs(U, coords, gp[2], 1/2, E, nu)
                            Ql = B*Chis[i][Ce,:]
                            Qr = B*Chis[j][Ce,:]
                            Ael[i,j] += gp[1]*(Ql'*A_2D*Qr)[1]*(dx/2)
                        end
                    end
                end
                return Ael
            end
        end

        A14 = zeros(Float64, 4, 4)
        for i in N:N:N*N
            Ce = vec(DOFs[i, :, :]')
            A_i = MacroA(1, U[Ce, :], i, elCoords[i, :, :], E, nu)
            A14 += A_i
        end

        A23 = zeros(Float64, 4, 4)
        for i in N*(N-1)+1:1:N*N
            Ce = vec(DOFs[i, :, :]')
            A_i = MacroA(2, U[Ce, :], i, elCoords[i, :, :], E, nu)
            A23 += A_i
        end

        return L/(L*L)*vcat(A14[1,:]', A23[2,:]', A23[3,:]', A14[4,:]')
    end
    return homogenizeP(), homogenizeA()
end