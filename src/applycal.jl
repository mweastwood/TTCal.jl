function applycal!(ms::Table,
                   gains::Array{Complex64,3};
                   force_imaging_columns::Bool = false)
    data = ms["DATA"]
    ant1,ant2 = ants(ms)
    applycal!(data,gains,ant1,ant2)
    if force_imaging_columns || Tables.checkColumnExists(ms,"CORRECTED_DATA")
        ms["CORRECTED_DATA"] = data
    else
        ms["DATA"] = data
    end
    data
end

function applycal!(data::Array{Complex64,3},
                   gains::Array{Complex64,3},
                   ant1::Vector{Int32},
                   ant2::Vector{Int32})
    Nbase = length(ant1)
    Nfreq = size(gains,3)
    for α = 1:Nbase, β = 1:Nfreq
        data[1,β,α] = (gains[ant1[α],1,β]*conj(gains[ant2[α],1,β]))\data[1,β,α]
        data[2,β,α] = (gains[ant1[α],1,β]*conj(gains[ant2[α],2,β]))\data[2,β,α] # this one could be
        data[3,β,α] = (gains[ant1[α],2,β]*conj(gains[ant2[α],1,β]))\data[3,β,α] # swapped with this one
        data[4,β,α] = (gains[ant1[α],2,β]*conj(gains[ant2[α],2,β]))\data[4,β,α]
    end
end

