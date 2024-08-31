using ArgCheck

# General configuration

timestep = 1. # 1 hour
SS = 36 # project for 1.5 days
# Noon to following midnight
mcdraws = 1 # 1 for deterministic

# Actions
PP = 11 # discretized power choices (including no-change action)

# States
#   vehicles_plugged: The number of plugged in cars
#   soc_plugged: The fraction of energy available for plugged-in cars
#   soc_driving: The fraction of energy available for driving cars
# Number of states for each dimension:
EE = 5 # 0 - 4 cars
FF = 11 # For both soc_plugged and soc_driving, 0 - 1

## Checks on configuration parameters

@argcheck mcdraws > 0
