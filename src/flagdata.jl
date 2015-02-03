################################################################################
# Public Interface

@doc """
Clear all of the flags in the measurement set.
""" ->
function clearflags!(ms::Table)
    N = numrows(ms)
    Nfreq = length(freq(ms))
    flags = zeros(Bool,4,Nfreq,N)
    ms["FLAG"] = flags
    flags
end

@doc """
Apply antenna flags to the measurement set.
""" ->
function flagdata!(ms::Table,flaggedantennas::Vector{Int})
    N = numrows(ms)
    Nfreq = length(freq(ms))
    ant1,ant2 = ants(ms)
    flags = ms["FLAG"]
    for α = 1:N
        if ant1[α] in flaggedantennas || ant2[α] in flaggedantennas
            flags[:,:,α] = true
        end
    end
    ms["FLAG"] = flags
    flags
end

