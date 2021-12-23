function createMesh(L::Float64, dx::Float64)
    N = Int(round(L/dx))
    Nc, Ni = 4, (N-1)^2
    Nm, Ns = 2*(N-1), 2*(N-1)
    Nel = N^2
    NDOF = 2*(N+1)*(N+1)
    nodes = zeros(Int64, N+1, N+1)
    
    # Corners
    nodes[end,1], nodes[end,end] = (1,2)
    nodes[1,end], nodes[1,1] = (3,4)
    # Internal
    nodes[2:N,2:N] = reshape(5 : Ni+4, (N-1, N-1))
    # Left
    nodes[2:end-1,1] = Ni+5 : Ni+N+3
    # Bottom
    nodes[end,2:end-1] = Ni+N+4 : Ni+2*(N+1)
    # Right
    nodes[2:end-1,end] = Ni+2*N+3 : Ni+3*N+1
    # Top
    nodes[1,2:end-1] = Ni+3*N+2 : Ni+4*N

    # Generate connectivity of each element based on their node numbers.
    NC = zeros(Int64, Nel, 4)
    for j in 1:N
        for i in 1:N
            NC[i+(j-1)*N, :] = [nodes[end-j+1,i], nodes[end-j+1,i+1], nodes[end-j,i+1], nodes[end-j,i]]
        end
    end
    # Create square mesh over the domain.
    x, y = (0:dx:L, 0:dx:L)
    XY = zeros(Float64, (N+1)^2, 2)
    for j in 1:N+1
        for i in 1:N+1
            XY[i+(j-1)*(N+1), :] = [x[i], y[j]]
        end
    end

    elCoords = zeros(Float64, Nel, 4, 2)
    for j in 1:N
        for i in 1:N
            k = i+(j-1)*(N+1)
            ind = [k, k+1, k+N+2, k+N+1]
            elCoords[i+(j-1)*N, :, :] = XY[ind,:]
        end
    end
    DOFs = zeros(Int64, Nel, 4, 2)
    for i in 1:Nel
        DOFs[i, :, :] = reshape(hcat(2*NC[i,:]' .- 1, 2*NC[i,:]'), (4,2))
    end

    # Dependency matrix
    rowInd = 1 : 2*Nm
    colInd = vcat(repeat([3,4], N-1), repeat([7,8], N-1))
    data = ones(Int64, 2*Nm)
    D = sparse(rowInd, colInd, data)
    Nci = Nc + Ni
    Nim = Ni + Nm
    Ncim = Nc + Ni + Nm
    Nims = Ni + Nm + Ns
    Ncims = Nc + Nims
    Ns = [Nc, Ni, Nm, Ns, Nel, NDOF, Nci, Nim, Ncim, Nims, Ncims]
    
   return N, NC, XY, elCoords, DOFs, D, Ns;
end