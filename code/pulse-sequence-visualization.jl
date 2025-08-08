function visualize_pulse_sequence!(plt, pulses::Vector{Tuple{Symbol, BigFloat}}; offset=0)
    ts = [-1,0.]
    zs = [0.,0]
    xs = [0.,0]
    for (d,t) in pulses
        push!(ts, ts[end])
        push!(ts, ts[end]+abs(t))
        s = copysign(0.9, t)
        if d == :z
            push!(zs, [s, s]...)
            push!(xs, [0,0]...)
        else
            push!(xs, [s,s]...)
            push!(zs, [0,0]...)
        end
    end
    push!(ts, ts[end])
    push!(ts, ts[end]+1)
    push!(xs, [0,0]...)
    push!(zs, [0,0]...)

    
    

    # plot sequence
    # z:
    Plots.plot!(plt, ts, zs .+ offset, c=:blue, lw=2, label="z")
    Plots.plot!(plt, ts, xs .+ 2 .+ offset, c=:red, lw=2, label="x")
    Plots.hline!(plt, [0 2] .+ offset, c=[:blue :red], ls=:dash, lw=1, label=false)
    #Plots.yaxis!(false)
    ytks = Plots.yticks(plt)[1]
    push!(ytks[1], ([0,2] .+ offset)...)
    push!(ytks[2], ["0", "0"]...)
    Plots.yticks!(plt, ytks[1], ytks[2])
    Plots.plot!(plt, xlims=extrema(ts), xlabel="time", size=(600, 200), legendposition=:outerright, margins=3*Plots.mm)
end

function visualize_pulse_sequence(pulses::Vector{Tuple{Symbol, BigFloat}})
    plt = Plots.plot()
    Plots.yticks!(plt, Real[], String[])
    visualize_pulse_sequence!(plt, pulses)
    return plt
end

function visualize_pulse_sequence!(pulses::Vector{Tuple{Symbol, BigFloat}}; offset=0)
    plt = Plots.plot!()
    visualize_pulse_sequence!(plt, pulses, offset=offset)
end