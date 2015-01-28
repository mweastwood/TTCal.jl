@doc """
The Butcher tableau for the 2nd-order Runge-Kutta method.
""" ->
const RK2_tableau = [1/2 0.0
                     0.0 1.0]

@doc """
The Butcher tableau for the 3rd-order Runge-Kutta method.
""" ->
const RK3_tableau = [1/2 0.0 0.0;
                    -1.0 2.0 0.0;
                     1/6 2/3 1/6]

@doc """
The Butcher tableau for the 4th-order Runge-Kutta method.
""" ->
const RK4_tableau = [1/2 0.0 0.0 0.0;
                     0.0 1/2 0.0 0.0;
                     0.0 0.0 1.0 0.0;
                     1/6 1/3 1/3 1/6]

@doc """
This singleton type is used to indicate which Runge-Kutta method
should be used. For example, `RK{4}` tells us to use the RK4 method.
""" ->
immutable RK{N}; end
const RK2 = RK{2}()
const RK3 = RK{3}()
const RK4 = RK{4}()

@doc """
This singleton type is used to pass the name of the step
function to be used by the Runge-Kutta method. This is preferrable
to directly passing the step function to `rkstep!` because
anonymous functions are slow. Instead, we can use a staged function
so that the compiler has knowledge of which function to call
at compile time (instead of runtime). Neat!

The step function must take the form:

    function func!(output,input,args)
        ...
    end

where `output` is overwritten with the next step to take given `input`.
""" ->
immutable RKInnerStep{func}; end

immutable RKWorkspace{T,N}
    x′::Vector{T}
    k::Vector{Vector{T}}
end
RKWorkspace{T}(x::Vector{T},N::Int) = RKWorkspace{T,N}(similar(x),[similar(x) for i = 1:N])

@doc """
Take a Runge-Kutta step.

* `x` is the initial point
* `x′` must be similar to `x` (same type and size)
* `k` is a vector of length `N` where each element is similar to `x` (arrays of the same type and size)
* `args` is simply passed as the third argument to `func!`

The output is stored in `x`, but all three variables are overwritten.
""" ->
stagedfunction rkstep!{func!,N}(x,x′,k,::RKInnerStep{func!},::RK{N},args...)
    tableau = symbol("RK$(N)_tableau")

    quote
        # Take the Runge-Kutta step
        $func!(k[1],x,args...)
        for row = 1:N-1
            for i = 1:length(x)
                x′[i] = x[i]
            end
            for col = 1:row
                $tableau[row,col] == 0.0 && continue
                for i = 1:length(x)
                    x′[i] += k[col][i]*$tableau[row,col]
                end
            end
            $func!(k[row+1],x′,args...)
        end

        # Output
        for col = 1:N
            for i = 1:length(x)
                x[i] += k[col][i] * $tableau[N,col]
            end
        end
        x
    end
end

function rkstep!{T,N}(x,workspace::RKWorkspace{T,N},_::RKInnerStep,args...)
    rkstep!(x,workspace.x′,workspace.k,_,RK{N}(),args...)
end

immutable StoppingCriteria
    maxiter::Int
    tolerance::Float64
end

