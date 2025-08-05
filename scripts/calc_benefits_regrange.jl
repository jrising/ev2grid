include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time
RR = 5 # number of possible regrange values
probfail_penalty = 10.
soc_plugged_1 = 0.5
soc_driving_1 = 0.5


benefits = run_optimized_regrange_simulation(dt0, SS, drive_starts_time, park_starts_time)
print("Optimized regrange benefits $(benefits)\n")

benefits_rot = thumbrule_regrange(drive_starts_time, park_starts_time)
print("Rule of thumb regrange benefits: $(benefits_rot)")