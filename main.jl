include("./helperFunctions.jl")
include("./Element.jl")
include("./Mesh.jl")
include("./NonlinearSolver.jl")
include("./LinearSolver.jl")
include("./Homogenizator.jl")

E, nu = 1.0, 0.3
L, dx = 2.0, 0.1

F = [[1.0 0.0 0.0]
     [0.0 1.5 0.0]
     [0.0 0.0 1.0]]
     
N, NC, XY, elCoords, DOFs, D, Ns = createMesh(L, dx)
println("Number of elements: ", N^2)

println("Solving nonlinear problem...")
U_nonlin = NonLinSolver(F, Ns, L, E, nu, 0.1)

mapI, mapJ = [1, 2, 1, 2], [1, 2, 2, 1]
Chis = [LinSolver(mapI[i], mapJ[i], U_nonlin) for i in 1:4]

println("Homogenization process...")
P, A = Homogenization(U_nonlin, Chis, N, E, nu)

println("Done!")
display(A)
readline()