function applycal(interferometer::Interferometer,
                  ms::Vector{Table},
                  gains)
    Nfreq = interferometer.Nfreq
    for i = 1:length(ms)
        applycal(interferometer,ms[i],gains[:,:,(i-1)*Nfreq+1:i*Nfreq])
    end
    nothing
end

function applycal(interferometer::Interferometer,
                  ms::Table,
                  gains)
    data = ms["DATA"]
    applycal!(data,gains)
    if Tables.checkColumnExists(ms,"CORRECTED_DATA")
        ms["CORRECTED_DATA"] = data
    else
        ms["DATA"] = data
    end
    nothing
end

function applycal(data::Array{Complex64,3},gains::Array{Complex64,3})
    corrected = copy(data)
    applycal!(corrected,gains)
    corrected
end

function applycal!(data::Array{Complex64,3},gains::Array{Complex64,3})
    Nant  = size(gains,1)
    Nfreq = size(gains,3)
    α = 1 # baseline counter
    for ant1 = 1:Nant, ant2 = ant1:Nant
        for β = 1:Nfreq
            data[1,β,α] = (gains[ant1,1,β]*conj(gains[ant2,1,β]))\data[1,β,α]
            data[2,β,α] = (gains[ant1,1,β]*conj(gains[ant2,2,β]))\data[2,β,α] # this one could be swapped
            data[3,β,α] = (gains[ant1,2,β]*conj(gains[ant2,1,β]))\data[3,β,α] # with this one
            data[4,β,α] = (gains[ant1,2,β]*conj(gains[ant2,2,β]))\data[4,β,α]
        end
        α += 1
    end
end

