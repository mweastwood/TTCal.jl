################################################################################
# Public Interface

function polcal(ms::Table,
                sources::Vector{Source},
                criteria::StoppingCriteria;
                force_imaging_columns::Bool = false,
                reference_antenna::Int = 1)
    frame = reference_frame(ms)
    u,v,w = uvw(ms)
    ν = freq(ms)
    ant1,ant2 = ants(ms)

    Nfreq = length(ν)
    Nant  = numrows(Table(ms[kw"ANTENNA"]))

    #data  = ms["DATA"]
    data  = ms["CORRECTED_DATA"] # We want data that is mostly right already
    model = genvis(frame,sources,u,v,w,ν)
    flags = ms["FLAG"]
    gains = ones(Complex64,2,2,Nant,Nfreq)

    if force_imaging_columns || Tables.checkColumnExists(ms,"MODEL_DATA")
        ms["MODEL_DATA"] = model
    end

    polcal!(gains,data,model,flags,ant1,ant2,criteria,reference_antenna)
    gains
end

################################################################################
# Internal Interface

function polcal!(gains::Array{Complex64,4},
                 data::Array{Complex64,3},
                 model::Array{Complex64,3},
                 flags::Array{Bool,3},
                 ant1::Vector{Int32},
                 ant2::Vector{Int32},
                 criteria::StoppingCriteria,
                 reference_antenna::Int)
    # Transpose the data and model arrays to create a better memory access pattern
    data  = permutedims(data, (1,3,2))
    model = permutedims(model,(1,3,2))
    flags = permutedims(flags,(1,3,2))

    Nfreq = size(gains,4)
    for β = 1:Nfreq
        polcal_onechannel!(slice(gains,:,:,:,β),
                           slice( data,:,:,β),
                           slice(model,:,:,β),
                           slice(flags,:,:,β),
                           ant1, ant2,
                           criteria,
                           reference_antenna)
    end
end

function polcal_onechannel!(gains, data, model, flags,
                            ant1::Vector{Int32},
                            ant2::Vector{Int32},
                            criteria::StoppingCriteria,
                            reference_antenna::Int)
    Nant = size(gains,3)

    square_data  = makesquare_polarized(data, ant1,ant2)
    square_model = makesquare_polarized(model,ant1,ant2)
    square_flags = makesquare_polarized(flags,ant1,ant2)

    # 2. Flag the auto-correlations.
    for i = 1:size(square_flags,1)
        square_flags[i,i] = true
    end

    # Start the gains off at something sensible.
    # This is an especially good approximation if we've
    # already applied a bandpass calibration.
    for i = 1:Nant
        gains[:,:,i] = [1 0;
                        0 1]
    end

    square_data  = square_data  .* !square_flags
    square_model = square_model .* !square_flags

    output = similar(gains)

    # Create rkstep! workspace variables.
    # Because "gains" is not a vector, it isn't serviced by the
    # current definition of RKWorkspace
    #workspace = RKWorkspace(gains,4)
    x′ = similar(gains)
    k  = [similar(gains) for i = 1:2]
    oldgains = similar(gains)

    # Iterate!
    iter = 0
    converged = false
    while !converged && iter < criteria.maxiter
        oldgains[:] = gains
        rkstep!(gains,x′,k,polcal_step!,RK2,square_data,square_model)
        if vecnorm(gains-oldgains)/vecnorm(oldgains) < criteria.tolerance
            converged = true
        end
        iter += 1
    end

    # Should fix the phase to something....
end

function makesquare_polarized(input,ant1::Vector{Int32},ant2::Vector{Int32})
    N = size(input,2)
    M = 2round(Integer,div(sqrt(1+8N)-1,2))
    output = zeros(eltype(input),M,M)
    for α = 1:N
        if ant1[α] == ant2[α]
            output[2ant1[α]-1,2ant1[α]-1] = input[1,α] # xx
            output[2ant1[α]-1,2ant1[α]-0] = input[2,α] # xy
            output[2ant1[α]-0,2ant1[α]-1] = input[3,α] # yx
            output[2ant1[α]-0,2ant1[α]-0] = input[4,α] # yy
        else
            output[2ant2[α]-1,2ant1[α]-1] = input[1,α] # x₁x₂
            output[2ant2[α]-0,2ant1[α]-1] = input[2,α] # x₁y₂
            output[2ant2[α]-1,2ant1[α]-0] = input[3,α] # y₁x₂
            output[2ant2[α]-0,2ant1[α]-0] = input[4,α] # y₁y₂

            output[2ant1[α]-1,2ant2[α]-1] = conj(input[1,α]) # x₂x₁
            output[2ant1[α]-1,2ant2[α]-0] = conj(input[2,α]) # y₂x₁
            output[2ant1[α]-0,2ant2[α]-1] = conj(input[3,α]) # x₂y₁
            output[2ant1[α]-0,2ant2[α]-0] = conj(input[4,α]) # y₂y₁
        end
    end
    output
end

function polcal_inner!(output,input,V,M)
    N = size(input,3)
    Z = similar(V)
    @inbounds for j = 1:N, i = 1:N
        # TODO: make this cleaner
        # The matrix multiplication is manually inlined to prevent the construction of
        # a temporary matrix. A cleaner way to express this would be nice.
        #Z[2i-1:2i,2j-1:2j] = sub(input,:,:,i)' * sub(M,2i-1:2i,2j-1:2j)
        Z[2i-1,2j-1] = conj(input[1,1,i]) * M[2i-1,2j-1] + conj(input[2,1,i]) * M[2i-0,2j-1]
        Z[2i-0,2j-1] = conj(input[1,2,i]) * M[2i-1,2j-1] + conj(input[2,2,i]) * M[2i-0,2j-1]
        Z[2i-1,2j-0] = conj(input[1,1,i]) * M[2i-1,2j-0] + conj(input[2,1,i]) * M[2i-0,2j-0]
        Z[2i-0,2j-0] = conj(input[1,2,i]) * M[2i-1,2j-0] + conj(input[2,2,i]) * M[2i-0,2j-0]
    end
    @inbounds for j = 1:N
        output[:,:,j] = sub(Z,:,2j-1:2j) \ sub(V,:,2j-1:2j) - sub(input,:,:,j)
    end
end
const polcal_step! = RKInnerStep{:polcal_inner!}()

