@doc """
Generate beam model to attenuate source model flux (to be used in genvis.jl).
As an initial attempt, uses the naive (1-cos(el)^1.6) beam, normalized to 1 at zenith.
""" ->
function beammodel(frame,phase_dir,l,m)
    dir   = lm2dir(phase_dir,l,m)
    az,el = dir2azel(frame,dir)

    att   = 1-cos(el)^1.6
    att
end


