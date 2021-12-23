function gaussPoints1D(n::Int64)
    # Gauss points
    if n == 2
        gpCoord = [-1/sqrt(3) 1/sqrt(3)]
        gpW = [1.0 1.0]
    elseif n == 3
        gpCoord = [-sqrt(3/5) 0 sqrt(3/5)]
        gpW = [5/9 8/9 5/9]
    elseif n == 4
        gpCoord = [-sqrt(3/7 + 2/7*sqrt(6/5)) -sqrt(3/7 - 2/7*sqrt(6/5)) sqrt(3/7 - 2/7*sqrt(6/5)) sqrt(3/7 + 2/7*sqrt(6/5))]
        gpW = [(18 - sqrt(30))/36 (18 + sqrt(30))/36 (18 + sqrt(30))/36 (18 - sqrt(30))/36]
    end
    return hcat(gpW', gpCoord')
end

function gaussPoints2D(n::Int64)
    # Gauss points
    coords, weights, x, w = zeros(n, n, 2), zeros(n, n), [], []
    if n == 2
        x = [-1 / sqrt(3), 1 / sqrt(3)]
        w = [1.0, 1.0]
    elseif n == 3
        x = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
        w = [5 / 9, 8 / 9, 5 / 9]
    elseif n == 4
        x = [-sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5)), -sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5)), sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5))]
        w = [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36]
    end
    for i in 1:n
        for j = 1:n
            coords[i, j, :] = [x[i], x[j]]
            weights[i, j] = w[i] * w[j]
        end
    end
    gpCoord = reshape(coords, n * n, 2)
    gpW = reshape(weights, n * n, 1)
    return hcat(gpW, gpCoord)
end