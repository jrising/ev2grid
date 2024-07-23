function periodstep(steps::Int)
    Hour(round(Int, timestep * steps)) + Minute((timestep * steps - round(Int, timestep * steps)) * 60)
end
