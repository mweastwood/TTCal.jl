function applycal(interferometer::Interferometer,
                  ms::Vector{MeasurementSet},
                  gains)
    Nfreq = interferometer.Nfreq
    for i = 1:length(ms)
        applycal(interferometer,ms[i],gains[:,:,(i-1)*Nfreq+1:i*Nfreq])
    end
    nothing
end

function applycal(interferometer::Interferometer,
                  ms::MeasurementSet,
                  gains)
    data = permutedims(getData(ms),(3,2,1))
    ant1 = getAntenna1(ms)
    ant2 = getAntenna2(ms)
    corrected = similar(data)

    Nfreq = interferometer.Nfreq
    Nbase = size(data,1)

    for β = 1:Nfreq, α = 1:Nbase
        corrected[α,β,1] = gains[ant1[α],1,β]\(conj(gains[ant2[α],1,β])\data[α,β,1])
        corrected[α,β,2] = gains[ant1[α],1,β]\(conj(gains[ant2[α],2,β])\data[α,β,2]) # this one could be swapped
        corrected[α,β,3] = gains[ant1[α],2,β]\(conj(gains[ant2[α],1,β])\data[α,β,3]) # with this one ?
        corrected[α,β,4] = gains[ant1[α],2,β]\(conj(gains[ant2[α],2,β])\data[α,β,4])
    end
    putCorrectedData!(ms,permutedims(corrected,(3,2,1)))
    #putData!(ms,permutedims(corrected,(3,2,1)))
    nothing
end

