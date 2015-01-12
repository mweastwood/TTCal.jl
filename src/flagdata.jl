################################################################################
# Public Interface

@doc """
Apply antenna flags to the measurement set.
""" ->
function flagdata!(ms::Table,interferometer::Interferometer)
    N = numrows(ms)
    # (the +1 converts to a 1-based indexing scheme)
    ant1 = ms["ANTENNA1"]+1
    ant2 = ms["ANTENNA2"]+1
    flags = zeros(Bool,4,channels(interferometer),N)
    for i = 1:N
        if ant1[i] in flaggedantennas(interferometer) || ant2[i] in flaggedantennas(interferometer)
            flags[:,:,i] = true
        end
    end
    ms["FLAG"] = flags
    flags
end

