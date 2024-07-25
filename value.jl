hourly_valuee = .1

## Value of having the energy at a given level
function value_energy(portion_below::Float64, enerfrac_avail::Float64, enerfrac_needed::Float64)
    if enerfrac_avail < enerfrac_needed
        return -Inf # Only happens if action made it so
    elseif enerfrac_avail > 1.
        return (1. - portion_below) - weight_portion_below * portion_below
    end

    return (1. - portion_below) * sqrt((enerfrac_avail - enerfrac_needed) / (enerfrac_max - enerfrac_needed)) - weight_portion_below * portion_below
end

value_energy(1., .3, .3)
value_energy(0., 0.5, 0.3)
value_energy(0., 0.5, 0.48)
value_energy(0., 0.5, 0.3)
value_energy(0., 0.6, 0.48)

function value_power_action(price::Float64, denerfrac::Float64)
    denergy = denerfrac * vehicle_capacity * vehicles # charging - discharging in terms of kWh

    if denergy > 0
        -price * denergy / efficiency # cost of energy
    else
        price * denergy # payment for energy
    end
end

function value_power_newstate(price::Float64, portion_below::Float64, enerfrac_toadd_below::Float64)
    denergy_below = enerfrac_toadd_below * vehicle_capacity * vehicles * portion_below
    return -price * denergy_below / efficiency
end
