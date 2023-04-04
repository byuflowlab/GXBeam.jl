function composite_pipe()

    E1 = 20.59e6
    E2 = 1.42e6
    nu12 = 0.42
    G12 = 0.87e6
    rho = 1.0
    pipemat = Material(E1, E2, E2, G12, G12, G12, nu12, nu12, nu12, rho)

    nx = 50
    nt = 24
    nr = 21
    nodes = Vector{Node{Float64}}(undef, (2*(nx-1) + 2*(nt-1))*nr)
    elements = Vector{MeshElement{Float64}}(undef, (2*(nx-1) + 2*(nt-1))*(nr-1))

    x1 = -1.0
    y1 = 0.3
    x2 = 1.0
    y2 = 0.5

    x = range(x1, x2, length=nx)
    y = range(y1, y2, length=nr)

    m = 1
    for i = 1:nx
        for j = 1:nr
            nodes[m] = Node(x[i], y[j])
            m += 1
        end
    end

    n = 1
    for i = 1:nx-1
        for j = 1:nr-1
            if j <= ((nr-1) ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = MeshElement([nr*(i-1)+j, nr*(i)+j, nr*(i)+j+1, nr*(i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I = nx-1

    r = y
    theta = reverse(range(-pi/2, pi/2, length=nt))
    for i = 2:nt
        for j = 1:nr
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]))
            m += 1
        end
    end

    for i = 1:nt-1
        for j = 1:nr-1
            if j <= ((nr-1) ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I += nt-1

    x1 = 1.0
    y1 = -0.3
    x2 = -1.0
    y2 = -0.5

    x = reverse(range(x2, x1, length=nx))
    y = reverse(range(y2, y1, length=nr))

    for i = 2:nx
        for j = 1:nr
            nodes[m] = Node(x[i], y[j])
            m += 1
        end
    end

    for i = 1:nx-1
        for j = 1:nr-1
            if j <= ((nr-1) ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I += nx-1

    theta = reverse(range(pi/2, 3*pi/2, length=nt))
    for i = 2:nt-1
        for j = 1:nr
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]))
            m += 1
        end
    end

    for i = 1:nt-1
        for j = 1:nr-1
            if i == nt-1
                Ip = -i
            else
                Ip = I
            end
            if j <= ((nr-1) ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(Ip+i)+j, nr*(Ip+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    return nodes, elements
end