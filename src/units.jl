const c = 2.99792e+8 * Meter/Second

# These functions are necessary because type inference
# currently fails such that [1Meter, 2Meter] != [1, 2]*Meter

@doc """
A convenience function for attaching units to an array.
""" ->
function addunits{T}(array::Array{T},unit::SIUnits.SIUnit)
    output = Array(quantity(T,unit),size(array))
    for i = 1:length(array)
        output[i] = array[i]*unit
    end
    output
end

@doc """
A convenience function for stripping units from an array.
""" ->
function stripunits(array::Array)
    output = Array(typeof(array[1].val),size(array))
    for i = 1:length(array)
        output[i] = array[i].val
    end
    output
end

