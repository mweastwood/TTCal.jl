"""
    multipeel!(visibilities, metadata, beam, sources)

Simultaneous removal of sources based on Smirnov & Tasse 2014.
"""
function multipeel!{T<:Source}(visibilities::Visibilities, meta::Metadata, beam::BeamModel,
                               sources::Vector{T}, N)
    frame = reference_frame(meta)
    sources = filter(sources) do source
        isabovehorizon(frame, source)
    end
    Nsource = length(sources)
    flags, visibilities_xx, visibilities_yy = reorder_visibilities(visibilities, meta)
    coherencies_xx, coherencies_yy = generate_coherencies(meta, beam, sources)
    gains_xx, gains_yy = generate_gains(meta, Nsource)
    gains_xx = multipeel!(visibilities_xx, coherencies_xx, gains_xx, N)
    # Output the gains
    output = [GainCalibration(Nant(meta), Nfreq(meta)) for source in sources]
    for β = 1:Nfreq(meta), ant = 1:Nant(meta), s = 1:length(sources)
        output[s].jones[ant,β] = DiagonalJonesMatrix(gains_xx[s, ant, β], gains_yy[s, ant, β])
    end
    for s = 1:length(sources)
        fixphase!(output[s], "1x")
    end
    output
end

function multipeel!(visibilities::Array{Complex128, 3}, coherencies::Array{Complex128, 4}, gains, N)
    #@code_warntype cohjones_step(gains, visibilities, coherencies)
    #for idx = 1:N
        #@time gains = cohjones_step(gains, visibilities, coherencies)
    #end
    #gains
    cohjones(gains, visibilities, coherencies, N)
end

function cohjones(gains, visibilities, coherencies, N)
    λ = 1.0 # Levenberg-Marquardt damping parameter
    ν = 2   # scaling factor for the damping parameter
    trial = similar(gains)
    best_χ² = cohjones_residual(gains, visibilities, coherencies)
    for idx = 1:N
        println("=====")
        @time next = cohjones_step(gains, visibilities, coherencies)

        attempts = 0
        successful = false
        if isodd(idx)
            λ = 0.0
        else
            λ = 1.0
        end
        #while !successful
            #attempts += 1
            @show λ

            A = λ/(1+λ)
            B = 1/(1+λ)
            for jdx in eachindex(gains)
                trial[jdx] = A*gains[jdx] + B*next[jdx]
            end

            χ² = cohjones_residual(trial, visibilities, coherencies)
            #if attempts
            #if χ² > best_χ²
            #    # unsuccessful step (try again with more damping)
            #    λ *= ν
            #elseif attempts > 1
            #    λ /= ν
            #    gains[:] = trial
            #    best_χ² = χ²
            #    successful = true
            #end
            @show χ²
            gains[:] = trial
        #end
        println("=====")
    end
    gains
end

function cohjones_residual(gains, visibilities, coherencies)
    output = 0.0
    Nant = size(visibilities, 1)
    Nfreq = size(visibilities, 2)
    Nsource = size(coherencies, 1)
    for ant2 = 1:Nant, β = 1:Nfreq, ant1 = 1:Nant
        δ = visibilities[ant1, β, ant2]
        for dir = 1:Nsource
            if ant1 == 1 && ant2 == 2 && β == 1
                @show gains[dir, ant1, β] gains[dir, ant2, β]
                @show coherencies[dir, ant1, β, ant2]
                @show visibilities[ant1, β, ant2]
            end
            δ -= gains[dir, ant1, β] * coherencies[dir, ant1, β, ant2] * conj(gains[dir, ant2, β])
        end
        output += abs2(δ)
    end
    output
end

function cohjones_step(input, visibilities, coherencies)
    Nant = size(visibilities, 1)
    Nfreq = size(visibilities, 2)
    Nsource = size(coherencies, 1)
    output = similar(input)
    JJ = zeros(Complex128, Nsource, Nsource)
    Jv = zeros(Complex128, Nsource)

    # As described by Smirnov & Tasse 2014, we need to compute two quantities:
    #  1. The Jacobian inner product with itself (JJ)
    #  2. The Jacobian inner product with the measured visibilities (Jv)
    # In both cases we do not need to explicitly compute the Jacobian itself. Instead we will use
    # the analytical expressions for the matrix elements.
    #
    # In doing this calculation we will assume that some (but not all!) of the off-diagonal elements
    # are zero. In particular we will include the off-diagonal terms corresponding to the
    # interactions between multiple directions. All other terms will be set to zero.

    for ant2 = 1:Nant, β = 1:Nfreq
        cohjones_setup(input, visibilities, coherencies, JJ, Jv, ant2, β)
        #@show ant2
        #println(JJ)
        #println(Jv)
        #if ant2 == 1 && β == 1
            #@show JJ Jv JJ\Jv
        #end
        cohjones_solve(output, JJ, Jv, ant2, β)
        #temp = cohjones_solve(JJ, Jv)
        #output[:, ant2, β] = temp
        #output[:, ant2, β] = bkfact(JJ)\Jv
    end
    output
end

function cohjones_setup(input, visibilities, coherencies, JJ, Jv, ant2, β)
    Nant = size(visibilities, 1)
    Nfreq = size(visibilities, 2)
    Nsource = size(coherencies, 1)
    JJ[:] = 0
    Jv[:] = 0
    @inbounds for ant1 = 1:Nant
        # Compute the contributions of the basline ant1-ant2

        # The Jacobian inner product with itself
        for dir2 = 1:Nsource
            y1 = input[dir2, ant1, β] * coherencies[dir2, ant1, β, ant2]
            for dir1 = dir2:Nsource
                y2 = input[dir1, ant1, β] * coherencies[dir1, ant1, β, ant2]
                JJ[dir1, dir2] += conj(y1) * y2
            end
            Jv[dir2] += y1 * conj(visibilities[ant1, β, ant2])
        end

        # The Jacobian inner product with the the measured visibilities
    end
    # This step might not be necessary if we cleverly specify which side of the matrix to use when
    # defining the Hermitian-ness
    @inbounds for dir2 = 1:Nsource, dir1 =  1:dir2-1
        JJ[dir1, dir2] = conj(JJ[dir2, dir1])
    end
    JJ, Jv
end

function cohjones_solve(output, JJ, Jv, ant2, β)
    factorization = bkfact!(JJ) # overwrites JJ
    #@time factorization = bkfact(JJ) # overwrites JJ
    A_ldiv_B!(view(output, :, ant2, β), factorization, Jv)
    #factorization \ Jv
end

function generate_coherencies(meta, beam, sources)
    initial_coherencies = [genvis(meta, beam, source) for source in sources]
    reorder_coherencies(initial_coherencies, meta)
end

function generate_gains(meta, Nsource)
    xx = ones(Complex128, Nsource, Nant(meta), Nfreq(meta))
    yy = ones(Complex128, Nsource, Nant(meta), Nfreq(meta))
    #sz = (Nsource, Nant(meta), Nfreq(meta))
    #xx = complex(randn(sz), randn(sz))
    #yy = complex(randn(sz), randn(sz))
    xx, yy
end

function reorder_visibilities(visibilities, meta)
    flags = zeros(Bool, Nant(meta), Nfreq(meta), Nant(meta))
    xx = zeros(Complex128, Nant(meta), Nfreq(meta), Nant(meta))
    yy = zeros(Complex128, Nant(meta), Nfreq(meta), Nant(meta))
    for β = 1:Nfreq(meta), α = 1:Nbase(meta)
        antenna1 = meta.baselines[α].antenna1
        antenna2 = meta.baselines[α].antenna2
        flags[antenna1, β, antenna2] = visibilities.flags[α, β]
        flags[antenna2, β, antenna1] = visibilities.flags[α, β]
        xx[antenna1, β, antenna2] = visibilities.data[α, β].xx
        xx[antenna2, β, antenna1] = visibilities.data[α, β].xx'
        yy[antenna1, β, antenna2] = visibilities.data[α, β].yy
        yy[antenna2, β, antenna1] = visibilities.data[α, β].yy'
    end
    flags, xx, yy
end

function reorder_coherencies(coherencies, meta)
    Nsource = length(coherencies)
    xx = zeros(Complex128, Nsource, Nant(meta), Nfreq(meta), Nant(meta))
    yy = zeros(Complex128, Nsource, Nant(meta), Nfreq(meta), Nant(meta))
    for s = 1:Nsource, β = 1:Nfreq(meta), α = 1:Nbase(meta)
        antenna1 = meta.baselines[α].antenna1
        antenna2 = meta.baselines[α].antenna2
        xx[s, antenna1, β, antenna2] = coherencies[s].data[α, β].xx
        xx[s, antenna2, β, antenna1] = coherencies[s].data[α, β].xx'
        yy[s, antenna1, β, antenna2] = coherencies[s].data[α, β].yy
        yy[s, antenna2, β, antenna1] = coherencies[s].data[α, β].yy'
    end
    xx, yy
end

