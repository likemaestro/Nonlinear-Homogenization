function Quad4(U::Matrix{Float64}, elID::Int64, coords::Matrix{Float64}, E::Float64, nu::Float64)
    elFint, elK = zeros(Float64, 8, 1), zeros(Float64, 8, 8)
    gps = gaussPoints2D(2)
    for gp in eachrow(gps)
        g, h = gp[2:end]
        B, detJ, PK1_1D, A_2D = constitutiveEqs(U, coords, g, h, E, nu)
        # Calculate the internal force matrix
        elFint += gp[1] * detJ * (B' * PK1_1D)
        # Calculate the element tangent stiffness matrix
        elK += gp[1] * detJ * (B' * A_2D * B)
    end
    return elFint, elK
end

function constitutiveEqs(U, coords, g, h, E, nu)
    Lambda, Mu = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))
    Nbar, B = zeros(Float64, 4, 2), zeros(Float64, 4, 8)
    PK1_1D, A_2D = zeros(Float64, 4, 1), zeros(Float64, 4, 4)
    mapI, mapJ = [1, 2, 1, 2], [1, 2, 2, 1]
    
    N_g = 1 / 4 * [-(1 - h) (1 - h) (1 + h) -(1 + h)]
    N_h = 1 / 4 * [-(1 - g) -(1 + g) (1 + g) (1 - g)]

    # Construct [Nbar]
    Nbar = [N_g' N_h']

    # Calculate Jacobian
    J = coords' * Nbar
    detJ, Jinv = det(J), inv(J)
    dNdX = Nbar * Jinv
    # Construct [B](4x8) matrix
    for i = 1:4
        B[1, 2*(i-1)+1] = dNdX[i, 1]
        B[2, 2*(i-1)+2] = dNdX[i, 2]
        B[3, 2*(i-1)+1] = dNdX[i, 2]
        B[4, 2*(i-1)+2] = dNdX[i, 1]
    end
    # Calculate deformation gradient
    F = reshape(U, 2, 4) * dNdX + I
    detF, Finv = det(F), inv(F)

    # First Piola-Kirchoff Stress Tensor
    PK1 = Mu * F - (Mu - Lambda * log(detF)) * Finv'

    # Nominal stiffness matrix
    A(i, j, k, l) = Mu * I[j, l] * I[i, k] + Lambda * Finv'[i, j] * Finv'[k, l] + (Mu - Lambda * log(detF)) * Finv'[i, l] * Finv'[k, j]

    # Rewrite PK1 and A in Voigt Notation
    for i = 1:4
        PK1_1D[i] = PK1[mapI[i], mapJ[i]]
        for j = 1:4
            A_2D[i, j] = A(mapI[i], mapJ[i], mapI[j], mapJ[j])
        end
    end
    return B, detJ, PK1_1D, A_2D
end