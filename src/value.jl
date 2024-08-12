hourly_valuee = .1


"""
Calculates the value of having energy at a certain level. The function considers the proportion of energy below
a threshold (`portion_below`), the available fraction of energy (`soc_avail`), and the needed fraction of
energy (`soc_needed`).

Returns:
- `-Inf` if the available energy is less than the needed energy, indicating an impossible state.
- A positive value considering the surplus energy if `soc_avail` is greater than 1,
  adjusted by `weight_portion_below`.
- A square-root scaled value when `soc_avail` is between `soc_needed` and `soc_max`,
  adjusted by the `weight_portion_below`.

Note: Ensure `soc_max` and `weight_portion_below` are defined in the scope where this function is used.
"""
function value_energy(portion_below::Float64, soc_avail::Float64, soc_needed::Float64)
    if soc_avail < soc_needed
        return -Inf # Only happens if action made it so
    elseif soc_avail > 1.
        return (1. - portion_below) - weight_portion_below * portion_below
    end

    return (1. - portion_below) * sqrt((soc_avail - soc_needed) / (soc_max - soc_needed)) - weight_portion_below * portion_below
end

"""
Calculates the value of an energy power action in terms of cost or payment. This is determined by the price
of energy (`price`) and the change in energy fraction (`dsoc`).

Returns:
- A negative cost of energy when `denergy` (charge discharge energy in kWh) is greater than 0, adjusted
  for charging efficiency.
- A payment received for energy when `denergy` is less than 0.

Note: Ensure `vehicle_capacity`, `vehicles`, and `efficiency` are defined in the scope where this function is used.
"""
function value_power_action(price::Float64, dsoc::Float64)
    denergy = dsoc * vehicle_capacity * vehicles # charging - discharging in terms of kWh

    if denergy > 0
        -price * denergy / efficiency # cost of energy
    else
        price * denergy # payment for energy
    end
end

"""
Calculates the cost associated with achieving a new energy state. This function evaluates the cost of adding
energy (`soc_toadd_below`) below a given portion (`portion_below`) at a specific energy price (`price`).

Returns:
- The cost of the additional energy needed below the current portion, adjusted for charging efficiency.

Note: Ensure `vehicle_capacity`, `vehicles`, and `efficiency` are defined in the scope where this function is used.
"""
function value_power_newstate(price::Float64, portion_below::Float64, soc_toadd_below::Float64)
    denergy_below = soc_toadd_below * vehicle_capacity * vehicles * portion_below
    return -price * denergy_below / efficiency
end
