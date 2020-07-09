module Tools

using LinearAlgebra
using Interpolations
using ..TimeMarching

export fluidbody_energy, fluidbody_momentum, curve_maxima, energy_harvested

"""
    fluidbody_momentum

Verify momentum conservation. Return change of body momentum through Δt + change
of fluid momentum - Lagrangian multipliers in the cdof of the first pinned joint
- mg. Currently only consider the revolute type as constrained type. Currently only
support 1d body moving in 2d space[x,y].
"""
function fluidbody_momentum(thist,bdhist,whist,fhist,λhist,bgs₀,Δx,xg,yg;gap=false,sample_rate=S1)

    # parameters
    Δt = thist[2] - thist[1]
    g = bdhist[1].sys.g
    nbody = bdhist[1].sys.nbody
    bgs = deepcopy(bgs₀)
    n = Int(size(thist,1)/sample_rate)
    cnt = 1

    # buffer
    trans_mom = zeros(n,2)
    t_sampled = zeros(n)
    λ6d = zeros(6)
    flu = zeros(n,2)
    flucur = zeros(n,2)
    flupre = zeros(n+1,2)
    kin = zeros(n,2)
    kincur = zeros(n,2)
    kinpre = zeros(n+1,2)
    lag = zeros(n,2)
    gra = zeros(n,2)

    for k = 1:sample_rate:size(thist,1)
        # time
        t_sampled[cnt] = thist[k]

        # --------- fluid momentum with a manually chosen origin ---------
        # # fluid momentum at current time
        # for (i,x) in enumerate(xg), (j,y) in enumerate(yg)
        #     flucur[cnt,1] += y*whist[k][i,j]
        #     flucur[cnt,2] -= x*whist[k][i,j]
        # end
        # flucur[cnt,:] .*= Δx
        # if k == 1
        #     flu[cnt,:] .= (flucur[cnt,:] - flupre[cnt,:])/Δt
        # else
        #     flu[cnt,:] .= (flucur[cnt,:] - flupre[cnt,:])/(sample_rate*Δt)
        # end
        # # record momentum into prev
        # flupre[cnt+1,:] .= flucur[cnt,:]

        # use fb from solutin directly
        ftmp = T₁ᵀ(bdhist[k],bgs,fhist[k],Δx; gap=gap)
        for i = 1:nbody
            ftmp2 = inv(bdhist[k].bs[i].Xb_to_i)'*ftmp[6i-5:6i]
            flu[cnt,1] -= ftmp2[4]
            flu[cnt,2] -= ftmp2[5]
        end

        # --------- body momentum ---------
        for i = 1:nbody
            # body momentum at current time
            tmp = inv(bdhist[k].bs[i].Xb_to_i)'*(bdhist[k].bs[i].inertia_b*bdhist[k].bs[i].v)
            kincur[cnt,:] .+= tmp[4:5]
        end
        if k == 1
            kin[cnt,:] .= (kincur[cnt,:] - kinpre[cnt,:])/Δt
        else
            kin[cnt,:] .= (kincur[cnt,:] - kinpre[cnt,:])/(sample_rate*Δt)
        end
        # record momentum into prev
        kinpre[cnt+1,:] .= kincur[cnt,:]

        # ---------  Lagrangian multipliers ---------
        if bdhist[1].js[1].joint_type == "revolute"
            λ6d[1:2] .= λhist[k][1:2]
            λ6d[4:6] .= λhist[k][3:5]
            lag[cnt,:] .+= (inv(bdhist[k].bs[1].Xb_to_i)'*λ6d)[4:5]
        end

        # --------- gravity ---------
        if g != zeros(Float64,3)
            for i = 1:nbody
                gra[cnt,2] -= bdhist[k].bs[i].mass*bdhist[k].sys.g[2]
            end
        end

        # --------- total ---------
        trans_mom[cnt,:] .= flu[cnt,:] .+ kin[cnt,:] .+ lag[cnt,:] .+ gra[cnt,:]
        cnt += 1
    end

    return t_sampled[2:end], trans_mom[2:end,:], flu[2:end,:], kin[2:end,:], lag[2:end,:], gra[2:end,:]
end


"""
    fluidbody_energy

Verify energy conservation. Return body kinetic energy + body gravitational potential
energy + spring potential energy + power extracted by damping*Δt - work done by fluid
on body*Δt - external work to enforce active motion*Δt.

For spring and gravitational potential energy, reference point is set as the initial
position. Currently only support 1d body moving in 2d space[x,y]
"""
function fluidbody_energy(thist,bdhist,whist,fhist,λhist,bgs₀,fsys;sample_rate=1,plane=[1,2])

    # parameters
    Δt = thist[2] - thist[1]
    Δx = fsys.Δx
    g = bdhist[1].sys.g
    bgs = deepcopy(bgs₀)
    nbody = bdhist[1].sys.nbody
    n = Int(size(thist,1)/sample_rate)
    cnt = 1

    # buffer
    energy = zeros(n)
    ulag = deepcopy(fhist[1])
    t_sampled = zeros(n)
    flu = zeros(n)
    dam = zeros(n)
    lag = zeros(n)
    kin = zeros(n)
    kincur = zeros(n)
    kinpre = zeros(n+1)
    spr = zeros(n)
    sprcur = zeros(n)
    sprpre = zeros(n+1)
    gra = zeros(n)
    gracur = zeros(n)
    grapre = zeros(n+1)

    for k = 1:sample_rate:size(thist,1)
        # time
        t_sampled[cnt] = thist[k]

        # --------- fluid energy ---------
        # get coordinates of Lagragian points
        for i = 1:length(bgs)
            b = bdhist[k].bs[bgs[i].bid]
            for j = 1:bgs[i].np
                q_temp = [zeros(Float64, 3); bgs[i].points[j]]
                q_temp = [zeros(Float64, 3); b.x_i] + b.Xb_to_i*q_temp
                bgs[i].q_i[j] = q_temp[4:6]
            end
        end
        fsys.X̃ = BodyGridToVectorData(bgs,"coord";plane=plane)
        # create B₂ operator based on current position of Lagrangian points
        regop_E = Regularize(fsys.X̃,Δx;weights=Δx^2,ddftype = Fields.Yang3)
        fsys.Emat = InterpolationMatrix(regop_E,fsys.Fq,fsys.Vb)
        ulag .= TimeMarching.B₂(whist[k],fsys)
        ulag.u .+= fsys.U∞[1]
        ulag.v .+= fsys.U∞[2]
        flu[cnt] = -Δx^2*dot(ulag,fhist[k])

        # --------- power extracted by damper ---------
        for i = 1:nbody
            vJtmp = bdhist[k].js[i].vJ
            for m = 1:bdhist[k].js[i].nudof
                dofid = bdhist[k].js[i].joint_dof[m].dof_id
                dam[cnt] += bdhist[k].js[i].joint_dof[m].damp*(vJtmp[dofid]^2)
            end
        end

        # --------- work by Lagrangian force to enforce active motion ---------
        atmp = 0 # temporary active dof
        for i = 1:nbody
            na = bdhist[k].js[i].na
            if na > 0
                udof_a = bdhist[k].js[i].udof_a
                anker = atmp
                for j = 1:na
                    anker = atmp + udof_a[j]
                    lag[cnt] += λhist[k][anker]*bdhist[k].js[i].vJ[udof_a[j]]
                end
            end
            atmp += 6-bdhist[k].js[i].np
        end

        # --------- body kinetic energy ---------
        for i = 1:bdhist[1].sys.nbody
            vtmp = bdhist[k].bs[i].v
            kincur[cnt] += 0.5*(bdhist[k].bs[i].inertia_b*vtmp)'*vtmp
        end
        if k == 1
            kin[cnt] = (kincur[cnt] - kinpre[cnt])/Δt
        else
            kin[cnt] = (kincur[cnt] - kinpre[cnt])/(sample_rate*Δt)
        end
        # record momentum into prev
        kinpre[cnt+1] = kincur[cnt]

        # --------- spring potential energy ---------
        for i = 1:nbody
            qJ0 = bdhist[1].js[i].qJ
            qJtmp = bdhist[k].js[i].qJ
            for m = 1:bdhist[k].js[i].nudof
                dofid = bdhist[k].js[i].joint_dof[m].dof_id
                # spring potential energy
                sprcur[cnt] += 0.5*bdhist[k].js[i].joint_dof[m].stiff*((qJtmp[dofid]-qJ0[dofid])^2)
            end
        end
        if k == 1
            spr[cnt] = (sprcur[cnt] - sprpre[cnt])/Δt
        else
            spr[cnt] = (sprcur[cnt] - sprpre[cnt])/(sample_rate*Δt)
        end
        # record momentum into prev
        sprpre[cnt+1] = sprcur[cnt]

        # --------- gravitational potential energy ---------
        if g != zeros(Float64,3)
            for i = 1:nbody
                y0 = 0.5*(bdhist[1].bs[i].verts_i[2,2] + bdhist[1].bs[i].verts_i[3,2])
                y = 0.5*(bdhist[k].bs[i].verts_i[2,2] + bdhist[k].bs[i].verts_i[3,2])
                gracur[cnt] += bdhist[k].bs[i].mass*abs(bdhist[k].sys.g[2])*(y-y0)
            end
        end
        if k == 1
            gra[cnt] = (gracur[cnt] - grapre[cnt])/Δt
        else
            gra[cnt] = (gracur[cnt] - grapre[cnt])/(sample_rate*Δt)
        end
        # record momentum into prev
        grapre[cnt+1] = gracur[cnt]

        # --------- total ---------
        energy[cnt] = flu[cnt] + dam[cnt] + lag[cnt] + gra[cnt] + kin[cnt] + spr[cnt]
        cnt += 1
    end

    return t_sampled[2:end], energy[2:end], flu[2:end], dam[2:end], lag[2:end],
            kin[2:end], spr[2:end], gra[2:end]
end

"""
    curve_maxima(thist,a1hist)

Given a periodic array, we want to identify if the curve is damping or oscillating
at a constant amplitude. This is done by finding out gradient of the curve and identify
local minima or maxima by change of gradient sign.

# Arguments

- `thist` : array of x-axis
- `a1hist` : array of y-axis
"""
function curve_maxima(thist,a1hist)
    # find out gradient of the curve using a1hist
    itp = interpolate(Array{Float64}(a1hist),BSpline(Linear()))
    g = [Interpolations.gradient(itp,i)[1] for i in 1:length(thist)]

    # identify local minimum or maximum by sign of gradient
    m = []
    for i in 1:length(thist)-1
        if sign(g[i])*sign(g[i+1]) == -1 push!(m,i) end
    end

    # print maximum and minimum
    output = a1hist[m]
    peak_time = thist[m]
    freq = []
    for i = 1:length(peak_time)-1
        push!(freq,0.5/(thist[m[i+1]]-thist[m[i]]))
    end
    println("peak amplitude 1: ",output[1:2:end],"\n")
    println("peak amplitude 2: ",output[2:2:end],"\n")
    println("frequency: ",freq)
    return output[1:2:end], output[2:2:end]

end

"""
    energy_harvested(t,q1,vJ1,q2,vJ2,c)

For a 2-linked body-joint system, calculate energy harvested through damper in one
period of oscillation, assuming both damper have the same damping coefficient c.

# Arguments

- `t` : history of time
- `q2` : history of joint 2's angle
- `vJ1` : history of joint 1's angular velocity
- `vJ2` : history of joint 2's angular velocity
- `c` : damping coefficient
"""
function energy_harvested(t,q1,q2,vJ1,vJ2,c)
    # find out gradient of the curve using q1
    itp1 = interpolate(Array{Float64}(q1),BSpline(Linear()))
    g1 = [Interpolations.gradient(itp1,i)[1] for i in 1:length(t)]
    # identify local minimum or maximum by sign of gradient
    m1 = []
    for i in 1:length(t)-1
        if sign(g1[i])*sign(g1[i+1]) == -1 push!(m1,i) end
    end
    period1 = t[m1[end]] - t[m1[end-2]]

    # find out gradient of the curve using q2
    itp2 = interpolate(Array{Float64}(q2),BSpline(Linear()))
    g2 = [Interpolations.gradient(itp2,i)[1] for i in 1:length(t)]
    # identify local minimum or maximum by sign of gradient
    m2 = []
    for i in 1:length(t)-1
        if sign(g2[i])*sign(g2[i+1]) == -1 push!(m2,i) end
    end

    # calculate harvested energy using time interval defined by m2
    period2 = t[m2[end]] - t[m2[end-2]]
    energy = 0.0
    for i = m2[end-2]:m2[end]
        energy += vJ1[i]^2 + vJ2[i]^2
    end
    energy  = energy*c/period2

    # output
    println("peak amplitude of joint 1: ",abs(q1[m1[end]]),"\n")
    println("peak amplitude of joint 2: ",abs(q2[m2[end]]),"\n")
    println("period of joint 1: ",period1,"\n")
    println("period of joint 2: ",period2,"\n")
    println("harvested energy per unit time: ",energy,"\n")
    println("harvested energy in one period: ",energy*period2,"\n")
    return
end

end # module end
