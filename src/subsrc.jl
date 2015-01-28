################################################################################
# Public Interface

function subsrc!{T<:AbstractSource}(ms::Table,sources::Vector{T})
    frame = reference_frame(ms)
    u,v,w = uvw(ms)
    ν     = freq(ms)
    data  = ms["CORRECTED_DATA"]

    subtracted = subsrc(frame,data,u,v,w,ν,sources)
    ms["CORRECTED_DATA"] = subtracted
    subtracted
end

################################################################################
# Internal Interface

function subsrc{T<:AbstractSource}(frame::ReferenceFrame,
                                   data::Array{Complex64,3},
                                   u::Vector{quantity(Float64,Meter)},
                                   v::Vector{quantity(Float64,Meter)},
                                   w::Vector{quantity(Float64,Meter)},
                                   ν::Vector{quantity(Float64,Hertz)},
                                   sources::Vector{T})
    model = genvis(frame,sources,u,v,w,ν)
    # Re-use the space allocated for the model visibilities
    # to store the model subtracted visibilities.
    for i = 1:length(model)
        model[i] = data[i]-model[i]
    end
    model
end

