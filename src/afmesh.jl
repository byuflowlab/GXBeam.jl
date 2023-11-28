# conveneince methods for initializing vectors and matrices
vector_ints(n) = Vector{Int64}(undef, n)

"""
    Layer(material, t, theta)

A layer (could be one ply or many plys of same material).
A layup is a vector of layers.

**Arguments**
- `material::Material`: material of layer
- `t::float`: thickness of layer
- `theta::float`: fiber orientation (rad)
"""
struct Layer{TF}
    material::Material{TF}
    t::TF
    theta::TF
end
Base.eltype(::Layer{TF}) where {TF} = TF
Base.eltype(::Type{Layer{TF}}) where {TF} = TF

Layer{TF}(l::Layer) where {TF} = Layer{TF}(l.material, l.t, l.theta)
Base.convert(::Type{Layer{TF}}, l::Layer) where {TF} = Layer{TF}(l)

# from https://discourse.julialang.org/t/findnearest-function/4143/4
function searchsortednearest(a, x)
    idx = searchsortedfirst(a, x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
       return idx
    else
       return idx-1
    end
end

"""
Modify number of layers based on a given maximum thickness
"""
function redistribute_thickness(segments, dt, nt)
    ns = length(segments)

    TF = eltype(eltype(eltype(segments)))
    newsegments = Vector{Vector{Layer{TF}}}(undef, ns)

    if !isnothing(nt)
        snt = sum.(nt)
        if !all(==(snt[1]), snt)
            error("number of elements in each subarray of nt must equal each other")
        end
        if any(vcat(nt...) .< 1)
            error("every element in nt must be greater than or equal to 1")
        end
    end

    for i = 1:ns
        # initialize new properties for this segment
        matvec = Material{TF}[]
        tvec = TF[]
        thetavec = TF[]

        for j = 1:length(segments[i])
            #extract layer
            layer = segments[i][j]
            # determine number of segments
            if isnothing(nt)
                nseg = round(Int64, layer.t/dt)
            else
                nseg = nt[i][j]
            end

            # decide whether to keep or divide up
            if nseg == 1 || nseg == 0  # keep this layer unchanged
                matvec = [matvec; layer.material]
                tvec = [tvec; layer.t]
                thetavec = [thetavec; layer.theta]
            else  # divide up existing layer into smaller pieces
                tvec = [tvec; fill(layer.t / nseg, nseg)]
                thetavec = [thetavec; fill(layer.theta, nseg)]  # copy over same theta and mat
                matvec = [matvec; fill(layer.material, nseg)]
            end
        end
        newsegments[i] = Layer.(matvec, tvec, thetavec)
    end

    return newsegments
end


"""
convert segments to all have the same number of layers for ease in meshing.
Overall definition remains consistent, just break up some thicker layers into multiple thinner layers (with same material and orientation properties)
"""
function preprocess_layers(segments, webs, dt=nothing, nt=nothing, wnt=nothing)

    TF = eltype(eltype(eltype(segments)))

    # number of segments
    ns = length(segments)

    # repartion thickneses if necessary so that thickness mesh is consistent
    if !isnothing(dt) || !isnothing(nt)
        segments = redistribute_thickness(segments, dt, nt)
    end
    if !isnothing(dt) || !isnothing(wnt)
        webs = redistribute_thickness(webs, dt, wnt)
    end

    # determine number of layers in each segment and thickness of each segment
    nl = vector_ints(ns)
    t = Vector{Vector{TF}}(undef, ns)
    for i = 1:ns
        nl[i] = length(segments[i])
        t[i] = Vector{TF}(undef, nl[i])
        for j = 1:nl[i]
            t[i][j] = segments[i][j].t
        end
    end


    # find segment with the most layers
    baseidx = argmax(nl)
    nlayers = nl[baseidx]

    # create new segments
    newsegments = Vector{Vector{Layer{TF}}}(undef, ns)

    for i = 1:ns
        # make all other segments have the same number of layers
        newsegments[i] = Vector{Layer{TF}}(undef, nlayers)

        # check if already has correct number of layers
        if length(segments[i]) == nlayers
            for j = 1:nlayers
                newsegments[i][j] = segments[i][j]
            end
            continue
        end

        # compute normalized cumulative thickness
        ti = cumsum(t[i])
        tb = cumsum(t[baseidx])
        tin = ti / ti[end]
        tbn = tb / tb[end]

        # map entries in tin to nerest value in tbn.  These thicknesses will remain unchanged
        for j = 1:length(tin)
            idx = searchsortednearest(tbn, tin[j])
            deleteat!(tbn, idx)
        end

        # append remaining values into original and resort
        new_tin = sort([tin; tbn])

        # convert back to thicknesses and find index of corresponding layer
        thickness = [new_tin[1]; diff(new_tin)] * ti[end]
        indices = Vector{Int64}(undef, nlayers)
        for j = 1:nlayers
            indices[j] = searchsortedfirst(tin, new_tin[j])
        end

        # create new segments using new thickness definitions
        for j = 1:nlayers
            seg = segments[i][indices[j]]
            newsegments[i][j] = Layer(seg.material, thickness[j], seg.theta)
        end

    end

    return newsegments, webs
end


"""
add new point xn into existing airfoil coordinates x, y.  Used to ensure breakpoints are in airfoil mesh.
assumes x is already sorted (and thus represents either upper or lower half)
"""
function insertpoint(x, y, xn)

    if indexin(xn, x) .== nothing  # element doesn't already exist in array
        idx = searchsortedfirst(x, xn)  # find where to insert
        yn = linear(x, y, xn)  # linear interpolate new y point

        # insert new point
        xnew = [x[1:idx-1]; xn; x[idx:end]]
        ynew = [y[1:idx-1]; yn; y[idx:end]]

        return xnew, ynew
    else
        return x, y
    end

end

"""
resample airfoil coordinates x, y with linear interpolation, making sure to exactly hit the points in vector xbreak (e.g., [0, 0.2, 0.4, 0.7, 1.0]).
ds defines the target size between coordinates (rounded to nearnest number of grid points).
assumes x is already sorted (and thus represents either upper or lower half or airfoil coordinates)
points in xbreak do not have to exist in the original x vector
alternatively you can define nseg, a vector of size: length(xbreak) - 1, defining how many points to use in each interval (e.g., [15, 20, 40, 30])
If nseg is provided then ds is ignored.
During optimization either keep mesh fixed, or if changes are needed, then nseg should be used rather than ds to avoid discrete changes.
"""
function resample(x, y, xbreak, ds, nseg=nothing)

    # initialize new airfoil coordinates
    xnew = [0.0]
    ynew = [0.0]

    # starting counter
    is = 1

    # iterate between each set of break points
    for i = 1:length(xbreak)-1

        # ending counter for where break exists (or will be inserted)
        ie = searchsortedfirst(x, xbreak[i+1])

        # extract start and end points of this interval
        xs = x[is]
        ys = y[is]
        xe = xbreak[i+1]  # use this rather than x[ie] in case xbreak does not already exist in x
        ye = linear(x[is:ie], y[is:ie], xe)

        # determine the number of points to use in this interval
        if isnothing(nseg)
            Deltas = sqrt((xe - xs)^2 + (ye - ys)^2)
            npts = round(Int64, Deltas/ds) + 1
        else
            npts = nseg[i] + 1
        end

        # create new segment and interpolate onto original airfoil coordinates
        if i === 1 || i === 2 || i === 4 #more resolution at LE
            xseg = halfcosine_spacing(xs, xe, npts)
            # xseg = range(xs, xe, length=npts)
            # xseg = halfsine_spacing(xs, xe, npts)
        elseif i === 3
            xseg = cosine_spacing(xs, xe, npts)
        elseif i === length(xbreak)-1 #more resolution at TE
            # xseg = halfcosine_spacing(xs, xe, npts)
            xseg = range(xs, xe, length=npts)
            # xseg = halfsine_spacing(xs, xe, npts)
        else #uniform spacing
            # xseg = halfcosine_spacing(xs, xe, npts)
            xseg = range(xs, xe, length=npts)
            # xseg = halfsine_spacing(xs, xe, npts)
        end
        yseg = linear(x[is:ie], y[is:ie], xseg)

        # add new segment to end of array
        xnew = [xnew; xseg[2:end]]
        ynew = [ynew; yseg[2:end]]

        # update starting index
        is = ie
    end

    return xnew, ynew
end


"""
airfoil coordinates xaf, yaf are assumed to start at trailing edge and traverse counterclockwise back to trailing edge.  (thus upper surface first)
coordinates are also assumed to be chord-normalized so that xaf is between 0 and 1.
xbreak is a vector defining locations that need to be inserted into airfoil coordinates (for upper and lower surfaces)
ds is a desired spacing for resampling (set ds=nothing to not resample)
Alternatively can use vector nseg to define number of points in each interval between xbreak (see resample())

Method separates airfoil into upper and lower surfaces, adds in break points, and optionally resamples.
"""
function parseairfoil(xaf, yaf, xbreak, ds, nseg=nothing)

    # ----- error checking -----
    if any(xaf .< 0.0) || any(xaf .> 1)
        error("airfoil coordinates must be chord normalized (0 <= x <= 1")
    end

    if xaf[1] != 1.0
        error("airfoil must start at trailing edge")
    end

    if xaf[end] != 1.0
        error("airfoil must end at trailing edge")
    end

    if yaf[3] < 0.0
        error("airfoil should traverse CCW (upper surface first)")
    end

    if xbreak[1] != 0.0
        error("xbreak must start at 0.0")
    end

    if xbreak[end] != 1.0
        error("xbreak must end at 1.0")
    end
    # -------------------------------

    # --- separate into upper and lower surfaces -----
    # find leading edge
    idx = argmin(xaf)
    le = xaf[idx]

    # check if leading edge corresponds to upper or lower surface
    if yaf[idx] > 0
        ustart = idx
        lstart = idx+1
    else
        ustart = idx-1
        lstart = idx
    end
    xu = xaf[ustart:-1:1]
    yu = yaf[ustart:-1:1]
    xl = xaf[lstart:end]
    yl = yaf[lstart:end]

    # force leading edge to start at 0 if not already (may slightly modify airfoil)
    if le != 0.0
        xu[1] = 0.0
        yu[1] = 0.0
        xl[1] = 0.0
        yl[1] = 0.0
    end

    # --- add in break points (xbreak) --------
    for i = 2:length(xbreak)-1  # omit start and end (which should be 0.0 and 1.0)
        xu, yu = insertpoint(xu, yu, xbreak[i])
    end
    for i = 2:length(xbreak)-1  # omit start and end (which should be 0.0 and 1.0)
        xl, yl = insertpoint(xl, yl, xbreak[i])
    end

    # ----- resample ---------
    if !isnothing(ds) || !isnothing(nseg)
        xu, yu = resample(xu, yu, xbreak, ds, nseg)
        xl, yl = resample(xl, yl, xbreak, ds, nseg)
    end

    return xu, yu, xl, yl
end



"""
determine tangential direction for each point on airfoil (as a normal vector with magnitudes tx and ty).
keyword signifies whether this is the upper or lower surface
"""
function tangential(x, y; upper=true)

    TF = promote_type(eltype(x), eltype(y))

    # initialize
    nx = length(x)
    tx = zeros(TF, nx)
    ty = zeros(TF, nx)

    # iterate through airfoil (except end points)
    for i = 2:nx-1
        # use central differencing to determine tangential direction
        dx = x[i+1] - x[i-1]
        dy = y[i+1] - y[i-1]
        ds = sqrt(dx^2 + dy^2)
        if upper
            tx[i] = dy/ds
            ty[i] = -dx/ds
        else
            tx[i] = -dy/ds
            ty[i] = dx/ds
        end
    end

    # tangent to leading always points straight back along chord line
    tx[1] = 1.0
    ty[1] = 0.0

    # tangent at last point alwyas point straight down along blunt T.E.
    tx[nx] = 0.0
    if upper
        ty[nx] = -1.0
    else
        ty[nx] = 1.0
    end

    return tx, ty
end

"""
determine tangential direction for each point on airfoil (as a normal vector with magnitudes tx and ty).  
keyword signifies whether this is the upper or lower surface
"""
function tangential_tyler(x, y, xbreak; upper=true)

    TF = promote_type(eltype(x), eltype(y))

    # initialize
    nx = length(x)
    tx = zeros(TF, nx)
    ty = zeros(TF, nx)

    x_central_LE = xbreak[2] * 1.5
    y_central_LE = xbreak[2] * 0.1894
    # y_central_LE = 0.0 #todo: find y centroid here
# println("$xbreak")
# println("x[1]: $(x[1])")
# println("x[end]: $(x[end])")
# println("y[1]: $(y[1])")
# println("y[end]: $(y[end])")
    # iterate through airfoil (except end points)
    for i = 1:nx-1
        idx = find_segment_idx(x[i], xbreak)

        if idx != 1 #previous method
            # use central differencing to determine tangential direction
            dx = x[i+1] - x[i-1]  
            dy = y[i+1] - y[i-1]
            ds = sqrt(dx^2 + dy^2)
            if upper
                tx[i] = dy/ds
                ty[i] = -dx/ds
            else
                tx[i] = -dy/ds
                ty[i] = dx/ds
            end
        else #new method
            dx = x_central_LE - x[i]
            dy = y_central_LE - y[i]
            ds = sqrt(dx^2 + dy^2)

            tx[i] = dx/ds
            ty[i] = dy/ds
        end
    end
    
    # tangent to leading always points straight back along chord line
    # tx[1] = 1.0
    # ty[1] = 0.0
    
    # tangent at last point alwyas point straight down along blunt T.E.
    tx[nx] = 0.0
    if upper
        ty[nx] = -1.0
    else
        ty[nx] = 1.0
    end 
    
    return tx, ty
end


function find_segment_idx(xi, xbreak)
    if xi == 0.0
        idx = 1
    else
        idx = searchsortedfirst(xbreak, xi) - 1
    end
    return idx
end


"""
determine curve for inner surface to find intersections at trailing edge and with webs
"""
function find_inner_surface(x, y, tx, ty, segments, xbreak)

    TF = promote_type(eltype(x), eltype(y), eltype(tx), eltype(ty))

    # initialize
    nx = length(x)

    xinner = Vector{TF}(undef, nx)
    yinner = Vector{TF}(undef, nx)
    for i = 1:nx
        # find corresponding segment
        idx = find_segment_idx(x[i], xbreak)

        # if idx != 1 #old method
            # pull out layers for this segment
            # local layers
            layers = segments[idx]
            nl = length(layers)  # number of layers

            # start at outer surface
            xinner[i] = x[i]
            yinner[i] = y[i]

            # add thicknesses of layers to reach inner surface
            for j = 1:nl
                xinner[i] += layers[j].t * tx[i]
                yinner[i] += layers[j].t * ty[i]
            end
        # else #if in 1st sector, try new method
        #     x0 = xbreak[2]
        #     y0 = 0.0
        #     xnew = -(x[i] - x0) #positive to the left of "origin"
        #     ynew = y[i]

        #     r = sqrt(xnew^2 + ynew^2)
        #     θ = atan(ynew,xnew)

        #     r2 = r - layers[j].t
        # end
    end

    return xinner, yinner
end



"""
determine where web intersects the inner surface curves.
creates a new discretization where points are added to line up with the vertical mesh in the web.
also returns the x vector index where the web mesh should start
"""
function web_intersections(xiu, yiu, xu, yu, txu, tyu, chord, webloc_i, web)

    TF = promote_type(eltype(xiu), eltype(yiu), eltype(xu), eltype(yu), eltype(txu), eltype(tyu))

    # compute total thickness of web
    webT = 0.0
    nl = length(web)
    for j = 1:nl
        webT += web[j].t
    end

    # center web on web loc
    xstart = webloc_i*chord - webT/2.0  # TODO: what if this is now out of the sector?  leave that up to user.

    # find x locations for the mesh discretizations through this web
    xiu_web = Vector{TF}(undef, nl+1)
    xiu_web[1] = xstart
    for j = 1:nl
        xiu_web[j+1] = xiu_web[j] + web[j].t
    end

    # compute corresponding points on airfoil surface
    yiu_web = Vector{TF}(undef, nl+1)
    xu_web = Vector{TF}(undef, nl+1)
    yu_web = Vector{TF}(undef, nl+1)
    for j = 1:nl+1
        # interpolate to find nondimensional distance
        idx = searchsortedlast(xiu, xiu_web[j])
        eta = (xiu_web[j] - xiu[idx]) / (xiu[idx+1] - xiu[idx])
        yiu_web[j] = yiu[idx] + eta * (yiu[idx+1] - yiu[idx])

        # project back onto airfoil surface
        xu_web[j] = xu[idx] + eta * (xu[idx+1] - xu[idx])
        yu_web[j] = yu[idx] + eta * (yu[idx+1] - yu[idx])
    end

    # compute the tangential direction for these new points
    deltax = xiu_web .- xu_web
    deltay = yiu_web .- yu_web
    mag = @. sqrt(deltax^2 + deltay^2)
    tx_web = deltax ./ mag
    ty_web = deltay ./ mag

    # find locations where to insert these new grid points (remove anything between them)
    idxs = searchsortedlast(xu, xu_web[1])
    idxe = searchsortedfirst(xu, xu_web[end])

    # create new grid points (note that other index values could change)
    newxiu = [xiu[1:idxs]; xiu_web; xiu[idxe:end]]
    newyiu = [yiu[1:idxs]; yiu_web; yiu[idxe:end]]
    newxu = [xu[1:idxs]; xu_web; xu[idxe:end]]
    newyu = [yu[1:idxs]; yu_web; yu[idxe:end]]
    newtxu = [txu[1:idxs]; tx_web; txu[idxe:end]];
    newtyu = [tyu[1:idxs]; ty_web; tyu[idxe:end]];

    return idxs+1, newxiu, newyiu, newxu, newyu, newtxu, newtyu
end


"""
given curves for inner surfaces, find where they intersect (if they intersect), for handling trailing-edge mesh
returns index of upper and lower surface where the trailing-edge-mesh should start (before intersection).
also returns the x, y location of intersection.
"""
function te_inner_intersection(xiu, yiu, xil, yil, xu, yu, xl, yl)

    # --- find crossing point -----
    # interpolate onto common x-grid
    yil2 = FLOWMath.linear(xil, yil, xiu)

    # subtract two curves
    ydiff = yiu - yil2

    # find first crossing on aft half of airfoil
    n = length(ydiff) ÷ 2  # integer division
    iu = findfirst(ydiff[n+1:end] .< 0.0)
    if isnothing(iu)  # no crossing
        return 0.0, 0.0, xu, yu, xl, yl
    end
    iu += n  # add back on number of points from fore half

    # find corresponding point on lower surface
    il = searchsortedfirst(xil, xiu[iu-1])
    # ----------------------------

    # 1 and 2 are endpoints of one line
    x1 = xiu[iu-1]; y1 = yiu[iu-1]
    x2 = xiu[iu]; y2 = yiu[iu]

    # 3 and 4 are endpoints of the other line
    x3 = xil[il-1]; y3 = yil[il-1]
    x4 = xil[il]; y4 = yil[il]

    # compute point of intersection between the two lines (https://dirask.com/posts/JavaScript-calculate-intersection-point-of-two-lines-for-given-4-points-VjvnAj)
    den = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
    px = (x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)
    py = (x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)

    # return point of intersection
    if den != 0
        px /= den
        py /= den

        # modify trailing edge mesh
        xu = xu[1:iu-1]; yu = yu[1:iu-1]
        xl = xl[1:il-1]; yl = yl[1:il-1]

        return px, py, xu, yu, xl, yl
    else  # no intersection
        return 0.0, 0.0, xu, yu, xl, yl
    end



end

"""
create nodes and elements for half (upper or lower) portion of airfoil
"""
function nodes_half(xu, yu, txu, tyu, xbreak, segments, chord, x_te, y_te)
    nl = length(segments[1])  # number of layers (same for all segments)

    TF = promote_type(eltype(xu), eltype(yu), eltype(txu), eltype(tyu), eltype(eltype(eltype(segments))), eltype(chord), eltype(x_te), eltype(y_te))
    # initialize
    nxu = length(xu)
    if x_te != 0.0
        nodes_te = nl + 1
        elements_te = nl
    else
        nodes_te = 0
        elements_te = 0
    end
    nodesu = Vector{Node{TF}}(undef, nxu * (nl + 1) + nodes_te)
    elementsu = Vector{MeshElement{TF}}(undef, (nxu - 1) * nl + elements_te)

    # create nodes
    n = 1
    # iterate across x
    for i = 1:nxu
        # start at airfoil surface
        x = xu[i]
        y = yu[i]
        nodesu[n] = Node(x, y)
        n += 1

        # find corresponding segment
        idx = find_segment_idx(xu[i], xbreak)
        layers = segments[idx]

        # check if at segment break
        segmentbreak = false
        if x == xbreak[idx+1] && xbreak[idx+1] != chord
            segmentbreak = true
        end

        for j = 1:nl
            # add thickness of each layer
            if !segmentbreak
                t = layers[j].t
            else
                t = 0.5 * (layers[j].t + segments[idx+1][j].t)  # average thickness at interface between two layers
            end

            x += t * txu[i]
            y += t * tyu[i]
            nodesu[n] = Node(x, y)
            n += 1
        end
    end

    # add trailing edge nodes
    if x_te != 0.0
        # pull out last row thickness and normalize it
        last_t = [layer.t for layer in segments[end]]
        norm_t = cumsum(last_t)
        norm_t /= norm_t[end]
        # distribute points according to that normalization starting at (chord, 0) ending at (x_te, y_te)
        nodesu[n] = Node(chord, zero(TF))
        n += 1
        for i = 1:length(norm_t)
            next_x = chord - (chord - x_te)*norm_t[i]
            next_y =  y_te*norm_t[i]
            nodesu[n] = Node(next_x, next_y)
            n += 1
        end
    else
        # modify last row of nodes to create angled transition from n-1 to n-nl
        # start
        nsx = nodesu[n-nl-2].x
        nsy = nodesu[n-nl-2].y
        # end
        nex = nodesu[n-nl-1].x
        ney = nodesu[n-nl-1].y

        for j = 1:nl+1
            njy = nodesu[n - (nl+1) - j].y  # use last row height (could really do either)

            frac = (njy - nsy)/(ney - nsy)
            nodesu[n - j] = Node(nsx + frac*(nex-nsx), njy)
        end

        # move second to last row back halfway to make room
        for j = 1:nl+1
            newx = 0.5 * (nodesu[n - (nl+1) - j].x + nodesu[n - 2*(nl+1) - j].x)
            newy = 0.5 * (nodesu[n - (nl+1) - j].y + nodesu[n - 2*(nl+1) - j].y)
            nodesu[n - nl-1 - j] = Node(newx, newy)
        end

        nxu -= 1  # hack for next part to create correct number of elements
    end

    # create elements
    e = 1
    for i = 1:nxu

        idx = find_segment_idx(xu[i], xbreak)
        layers = segments[idx]

        for j = 1:nl
            elementsu[e] = MeshElement([(nl+1)*(i-1) + (j+1), (nl+1)*i + (j+1), (nl+1)*i + j, (nl+1)*(i-1) + j], layers[j].material, layers[j].theta)
            e += 1
        end
    end

    return nodesu, elementsu
end

"""
given nodes/elements for upper surface and lower surface separately,
combine into one set making sure to reuse the common nodes that occur at the LE/TE
"""
function combine_halfs(nodesu, elementsu, nodesl, elementsl, nlayers, x_te)

    TN = promote_type(eltype(eltype(nodesu)), eltype(eltype(nodesl)))
    TE = promote_type(eltype(eltype(elementsu)), eltype(eltype(elementsl)))

    nt = 1 + nlayers  # number of points across thickness
    nnu = length(nodesu)
    nnl = length(nodesl)
    neu = length(elementsu)
    nel = length(elementsl)
    if x_te != 0.0
        nn = nnu + nnl - 2*nt  # number of nodes
        ne = neu + nel  # number of elements
    else
        nn = nnu + nnl - nt  # no shared t.e.
        ne = neu + nel + nlayers
    end

    nodes = Vector{Node{TN}}(undef, nn)
    elements = Vector{MeshElement{TE}}(undef, ne)

    # copy over upper nodes and elements unchanged
    nodes[1:nnu] .= nodesu
    elements[1:neu] .= elementsu

    # for lower surface we share leading and traling edges so don't copy over those nodes
    for i = nnu+1:nn
        j = i - nnu + nt  # starts at 1 + nt
        nodes[i] = Node(nodesl[j].x, nodesl[j].y)
    end

    # we retain the same number of elements, but the node numbers have changed on the lower surface
    # first nt-1 elements use the node numbers from the upper surface leading edge
    for i = neu+1:neu+nt-1
        j = i - neu  # starts at 1
        nodenum = [nnu+j+1; j+1; j; nnu+j]
        elements[i] = MeshElement(nodenum, elementsl[j].material, elementsl[j].theta)
    end

    # middle section need to reorder nodes
    if x_te != 0.0
        lastidx = neu + nel - (nt-1)
    else
        lastidx = neu + nel  # no shared t.e.
    end
    for i = neu+nt:lastidx
        j = i - neu  # element index (starts at nt)
        oldnodenum = elementsl[j].nodenum
        oldnodenum .+= (nnu - nt)  # increase node number by upper surface (minus the ones reused on leading edge)
        nodenum = [oldnodenum[2]; oldnodenum[1]; oldnodenum[4]; oldnodenum[3]]  # reorder the old nodenumbers since lower surface is flipped (sttart at bottom left, ccw)
        elements[i] = MeshElement(nodenum, elementsl[j].material, elementsl[j].theta)
    end

    # last nt-1 elements use the node numbers from the upper surface trailing edge
    if x_te != 0.0
        for i = neu+nel-(nt-1)+1:neu+nel
            j = i - neu  # element index
            k = i - (neu+nel-(nt-1))  # starts at 1
            oldnodenum = elementsl[j].nodenum
            oldnodenum .+= (nnu - nt)  # increase node number by upper surface (minus the ones reused on leading edge)
            nodenum = [nnu-nt+k+1; oldnodenum[1]; oldnodenum[4]; nnu-nt+k]
            elements[i] = MeshElement(nodenum, elementsl[j].material, elementsl[j].theta)
        end
    else
        # add new elements to close trailing edge.
        for i = 1:nlayers
            nodenum = [nnl-i+1+(nnu-nt); nnl-i+(nnu-nt); nnu-i; nnu-i+1]
            elements[neu+nel+i] = MeshElement(nodenum, elementsu[end-i+1].material, elementsu[end-i+1].theta)
        end
    end

    return nodes, elements
end

"""
add the nodes and elements for the webs.  given the x locations (by idx) where the webs start
and the number of grid points in the webs.
"""
function addwebs(idx_webu, idx_webl, nx_web, nodes, elements, webs, nnu, nl, ne_web=4)
    nt = 1 + nl  # number of points across thickness

    TN = eltype(eltype(nodes))
    TE = promote_type(eltype(eltype(elements)), eltype(eltype(eltype(webs))))

    # find nodes numbers for start of webs
    idx_webu *= nt  # there are nt points per nx index
    idx_webl *= nt
    idx_webl .+= (nnu - nt)  # lower surface node count is offset by this many nodes

    # initialize sizes of nodes and elements for web
    web_nodes = Vector{Node{TN}}(undef, (ne_web-1)*sum(nx_web))
    web_elements = Vector{MeshElement{TE}}(undef, ne_web*sum(nx_web .- 1))
    nn = length(nodes)
    n_web = 1
    e_web = 1

    # create nodes in webs
    for i = 1:length(nx_web)  # for each web
        for j = 1:nx_web[i]  # for each x direction in this web
            for k = 1:ne_web-1  # for each vertical direction in this x location
                nodeu = nodes[idx_webu[i] + (j-1)*nt]
                nodel = nodes[idx_webl[i] + (j-1)*nt]
                x = nodeu.x  # should be same top and bottom
                y = nodeu.y - k/ne_web*(nodeu.y - nodel.y)
                web_nodes[n_web] = Node(x, y)
                n_web += 1
            end
        end
    end

    # create elements
    for i = 1:length(nx_web)  # for each web
        start = nn + (i-1)*nx_web[i]*(ne_web-1)
        for j = 1:nx_web[i]-1  # for each x direction in this web
            for k = 1:ne_web  # for each vertical direction in this x location
                if k == 1
                    nodenum = [start + (ne_web-1)*(j-1) + k; start + (ne_web-1)*(j) + k; idx_webu[i] + (j)*nt; idx_webu[i] + (j-1)*nt]
                elseif k == ne_web
                    nodenum = [idx_webl[i] + (j-1)*nt; idx_webl[i] + (j)*nt; start + (ne_web-1)*(j) + k-1; start + (ne_web-1)*(j-1) + k-1]
                else
                    start = nn + (i-1)*nx_web[i]*(ne_web-1)
                    nodenum = [start + (ne_web-1)*(j-1) + k; start + (ne_web-1)*(j) + k; start + (ne_web-1)*(j) + k-1; start + (ne_web-1)*(j-1) + k-1]
                end
                web_elements[e_web] = MeshElement(nodenum, webs[i][j].material, webs[i][j].theta)
                e_web += 1
            end
        end
    end

    # concatenate
    nodes = [nodes; web_nodes]
    elements = [elements; web_elements]

    return nodes, elements
end


"""
    afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs; ds=nothing, dt=nothing, ns=nothing, nt=nothing, wns=4, wnt=nothing)

Create structural mesh for airfoil.  The airfoil coordinates define the meshing density tangential to the airfoil.
Whereas the number of layers defines the resolution normal to the airfoil. All segments are meshed with the same resolution
in the normal direction, using the number of grid points as defined by segment with the most layers.

**Arguments**
- `xaf, yaf::Vector{float}`: points defining airfoil start at trailing edge and traverse counterclockwise back to trailing edge (can be blunt or sharp T.E.)
- `chord::float`: chord length
- `twist::float`: twist angle (rad)
- `paxis::float`: pitch axis (normalized by chord). e.g., 0.25 means twist about quarter chord.
- `xbreak::Vector{float}`: x-locations, normalized by chord, defining break points between segments. must start at 0 and end at 1. e.g., [0, 0.2, 0.4, 0.7, 1.0] defines 4 segments.
- `webloc::Vector{float}`: x-locations, normalized by chord, defining web centers (and length of vector is number of webs).  e.g., [0.25, 0.55] means there is a web at 25% chord and a second web at 55% chord.
- `segments::Vector{Vector{Layer}}`: A Layer defines a ply (or multiple plys with the same material/orientation).  At a given x location the ply stack (segment) is defined a vector of layers starting from outside surface towards inside. Segments then is a vector of these segments that defines properties between segments as defined by xbreak.
- `webs::Vector{Vector{Layer}}`: same structure as segments, except each inner vector is from left to right (although this is usually symmetric), and each outer vector is for a separate web
- `ds::float`: if provided, airfoil spacing will be resampling with approximately this spacing, normalized by chord.  e.g., 0.01 will have points on the airfoil roughly 1% chord apart.
- `dt::float`: if provided, thickness will be resampled with this maximum mesh size (thickness is absolute). Note that the total number of cells remains constant along airfoil, so most thicknesses will be much less.  e.g., 0.01 will target a maximum mesh thickness of 0.01 (absolute).
- `ns::vector{int}`: if provided, rather than use a targert size ds, we specify the number of cells to use in each segment.  This is desirable for gradient-based optimization, if airfoil coordinates are changed, so that during resizing operations the  mesh stretch/shrinks rather than experiencing discrete jumps.  For example, ns=[15, 20, 40, 30] would use 15 elements between xbreak[1] and xbreak[2] and so on.
- `nt::vector{vector{int}}`: if provided, defines how many elements to use across tangential direction.  Again, prefered over dt for gradient-based optimization, if the thicknesses are changed during optimization.  each entry defines how many cells to put in that layer following order of original layup.  for example, nt=[[1, 2, 1], [1, 3]] would use 1 element, 2 elements (subdivide), then 1 elements over first sector, and so on.
- `wns::int`: discretization level for number of elements vertically along web.
- `wnt::vector{vector{int}}`: same definition as nt but for the webs

**Returns**
- `nodes::Vector{Node{Float64}}`: nodes for this mesh
- `elements::Vector{MeshElement{Float64}}`: elements for this mesh
"""
function afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs; ds=nothing, dt=nothing, ns=nothing, nt=nothing, wns=4, wnt=nothing)

    # -------------- preprocessing -----------------
    # preprocess the segments so all have same number of layers
    segments, webs = preprocess_layers(segments, webs, dt, nt, wnt)

    # separate into upper and lower surfaces
    xu, yu, xl, yl = parseairfoil(xaf, yaf, xbreak, ds, ns)

    # unnormalize
    xu *= chord; yu *= chord
    xl *= chord; yl *= chord
    xbreak *= chord

    # compute tangential directions
    txu, tyu = tangential_tyler(xu, yu, xbreak, upper=true)
    txl, tyl = tangential_tyler(xl, yl, xbreak, upper=false)
    # txu, tyu = tangential(xu, yu, upper=true)
    # txl, tyl = tangential(xl, yl, upper=false)
    
    # find inner surface curves
    xiu, yiu = find_inner_surface(xu, yu, txu, tyu, segments, xbreak)
    xil, yil = find_inner_surface(xl, yl, txl, tyl, segments, xbreak)

    # add webs. note that doing so changes the mesh so tangential directions and inner surface must be recomputed
    # webloc[i] must be in increasing order so that the previous indices in idx_web remain correct as points are added behind them.
    nw = length(webs)
    idx_webu = vector_ints(nw)
    idx_webl = vector_ints(nw)
    nx_web = vector_ints(nw)
    for i = 1:nw
        idx_webu[i], xiu, yiu, xu, yu, txu, tyu = web_intersections(xiu, yiu, xu, yu, txu, tyu, chord, webloc[i], webs[i])
        idx_webl[i], xil, yil, xl, yl, txl, tyl = web_intersections(xil, yil, xl, yl, txl, tyl, chord, webloc[i], webs[i])
        nx_web[i] = length(webs[i]) + 1
    end

    # determine intersection point for trailing edge.  (note must be done at end)
    x_te, y_te, xu, yu, xl, yl = te_inner_intersection(xiu, yiu, xil, yil, xu, yu, xl, yl)
    # -----------------------------------------------------------------

    # ------------------ build mesh --------------------
    nodesu, elementsu = nodes_half(xu, yu, txu, tyu, xbreak, segments, chord, x_te, y_te)
    nodesl, elementsl = nodes_half(xl, yl, txl, tyl, xbreak, segments, chord, x_te, y_te)

    nlayer = length(segments[1])
    nodes, elements = combine_halfs(nodesu, elementsu, nodesl, elementsl, nlayer, x_te)

    if nw > 0 # only add webs if there are webs defined
        nodes, elements = addwebs(idx_webu, idx_webl, nx_web, nodes, elements, webs, length(nodesu), nlayer, wns)
    end
    # -----------------------------------

    # ------ rotate with twist -------
    c = cos(twist)
    s = sin(twist)
    xc = paxis * chord
    for i = 1:length(nodes)
        x = nodes[i].x
        y = nodes[i].y
        nodes[i] = Node((x - xc)*c + y*s + xc, -(x - xc)*s + y*c)  # added xc back to keep origin at original location
    end
    # -----------------------------------

    return nodes, elements
end


function halfcosine_spacing(x1, x2, npoints)

    θs = range(0, pi/2, length=npoints)

    return (x2-x1) * (1 .- cos.(θs)) .+ x1
end

halfcosine_spacing(x::AbstractVector) = halfcosine_spacing(x[1], x[end], length(x))

function cosine_spacing(x1, x2, npoints)

    θs = range(0, pi, length=npoints)

    return (x2-x1) / 2 * (1 .- cos.(θs)) .+ x1
end

cosine_spacing(x::AbstractVector) = cosine_spacing(x[1], x[end], length(x))

function halfsine_spacing(x1, x2, npoints)

    θs = range(0, pi/2, length=npoints)

    return (x2-x1) * (sin.(θs)) .+ x1
end

halfsine_spacing(x::AbstractVector) = halfsine_spacing(x[1], x[end], length(x))

function sine_spacing(x1, x2, npoints)

    θs = range(0, pi, length=npoints)

    return (x2-x1) / 2 * (sin.(θs)) .+ x1
end

sine_spacing(x::AbstractVector) = sine_spacing(x[1], x[end], length(x))
