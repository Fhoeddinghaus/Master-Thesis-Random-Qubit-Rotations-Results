using ProgressMeter

PI = BigFloat(π)

# Metropolis step for a given function f_min(ωs, pulses)
function metropolis_min_step!(
        f_min::Function,       # Function(ωs, pulses) to minimize
        ωs::Vector{BigFloat}, 
        pulses::Vector{Tuple{Symbol, BigFloat}}, 
        fv_last::BigFloat,     # last value of f_min
        T::BigFloat, 
        stepsize::BigFloat=big(0.1),
        idx=false
    )
    # choose random step in sequence
    if idx == false
        idx = 1:length(pulses)
    end
    i = rand(idx) 

    pulses_new = copy(pulses)

    dir, t = pulses_new[i]
    min_ω, max_ω = extrema(ωs)
    
    t = (2*rand(BigFloat)-1) * stepsize + t
    if dir == :z
        # multiple rotations could be possible, 
        # but maybe limit it by the minimal velocity to 2-3 rotations?
        ϕ = (min_ω * t) % (2*2PI)
        t = ϕ/min_ω
    else
        # rot around x is limited
        t = t % (2PI)
    end
    pulses_new[i] = (dir, t)
    fv_new = f_min(ωs, pulses_new)
    delta = fv_new - fv_last

    if rand(BigFloat) < exp(-delta/T)
        pulses[i] = (dir, t)
        return fv_new
    end
    return fv_last
end

# Simulated Annealing for minimization of a function f_min(ωs, pulses)
function simulate_annealing_min(
        f_min::Function,       # Function(ωs, pulses) to minimize
        Ts::Vector{BigFloat}, 
        ωs::Vector{BigFloat}, 
        pulses_start::Vector{Tuple{Symbol, BigFloat}}; 
        sweepsize::Int64=10, 
        stepsize::Function=((i,T, fv_last) -> 0.1),
        idx=false,
        max_order_diff=3,
        max_order=9,
        p_desc="Simulating...",
        pulse_history=false
    )
    
    Fs = zeros(BigFloat, length(Ts) + 1) 
    fval = f_min(ωs, pulses_start) 
    Fs[1] = fval
    step = zeros(BigFloat, length(Ts)) 
    pulses = copy(pulses_start)
    
    if pulse_history isa Vector
        push!(pulse_history, copy(pulses)) 
    end

    last_index = length(Ts) + 1
    
    p = Progress(length(Ts), p_desc, showspeed=true)
    for (i,T) in enumerate(Ts) 
        step[i] = stepsize(i,T, fval)
        for _ in 1:sweepsize
            fval = metropolis_min_step!(f_min, ωs, pulses, fval, T, step[i], idx)
        end
        Fs[i+1] = fval
        
        if pulse_history isa Vector
            push!(pulse_history, copy(pulses)) 
        end
        
        next!(p, 
            showvalues=[
                (:iter, i),
                (:T, T),
                (:fval, fval),
                (:order, Float64(round(log10(fval), digits=2)))
            ]
        )
        # compare orders of T and fval
        #if log10(T)/log10(fval) >= max_order_diff
        if (log10(fval) - log10(T) >= max_order_diff) || (log10(fval) <= -max_order)
            println("Stopped @ (i,T, order) = ($i,$T, $(Float64(round(log10(fval), digits=2)))) with val = $fval")
            finish!(p,
                showvalues=[
                    (:iter, i),
                    (:T, T),
                    (:fval, fval),
                    (:order, Float64(round(log10(fval), digits=2)))
                ]
            )
            last_index = i + 1
            break
        end
    end

    return Fs[1:last_index], pulses, step
end

# Simulated Annealing for minimization of a function f_min(ωs, pulses) and a variable temperature T(i,T_true,fv_last)
function simulate_annealing_min_variable_T(
        f_min::Function,       # Function(ωs, pulses) to minimize
        Ts::Vector{BigFloat},  # controls maximal number of iterations
        ωs::Vector{BigFloat}, 
        pulses_start::Vector{Tuple{Symbol, BigFloat}}; 
        sweepsize::Int64=10, 
        stepsize::Function=((i,T, fv_last) -> 0.1),
        Tfunc::Function=((i,T,fv_last) -> T),
        idx=false,
        max_order_diff=3,
        max_order=9,
        p_desc="Simulating...",
        pulse_history=false
    )
    
    Fs = zeros(BigFloat, length(Ts) + 1) 
    fval = f_min(ωs, pulses_start) 
    Fs[1] = fval
    step = zeros(BigFloat, length(Ts)) 
    pulses = copy(pulses_start)
    
    if pulse_history isa Vector
        push!(pulse_history, copy(pulses)) 
    end

    last_index = length(Ts) + 1
    
    p = Progress(length(Ts), p_desc, showspeed=true)
    for (i,T) in enumerate(Ts) 
        T_ = T
        T = Tfunc(i, T_, fval)
        step[i] = stepsize(i,T_, fval)
        for _ in 1:sweepsize
            fval = metropolis_min_step!(f_min, ωs, pulses, fval, T, step[i], idx)
        end
        Fs[i+1] = fval
        
        if pulse_history isa Vector
            push!(pulse_history, copy(pulses)) 
        end
        
        next!(p, 
            showvalues=[
                (:iter, i),
                (:T, T),
                (:T_, T_),
                (:fval, fval),
                (:order, Float64(round(log10(fval), digits=2)))
            ]
        )
        # compare orders of T and fval
        #if log10(T)/log10(fval) >= max_order_diff
        if (log10(fval) - log10(T) >= max_order_diff) || (log10(fval) <= -max_order)
            println("Stopped @ (i,T,T_, order) = ($i,$T, $T_, $(Float64(round(log10(fval), digits=2)))) with val = $fval")
            finish!(p,
                showvalues=[
                    (:iter, i),
                    (:T, T),
                    (:fval, fval),
                    (:order, Float64(round(log10(fval), digits=2)))
                ]
            )
            last_index = i + 1
            break
        end
    end

    return Fs[1:last_index], pulses, step
end