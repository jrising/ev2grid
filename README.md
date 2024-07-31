# ev2grid
Vehicle-to-grid optimization

This algorithm operates at the aggregator scale, and is only used for
day-ahead planning:
 - Scheduled EV actions are aggregated
 - Instant charge is modeled as removing an EV from the pool
 
The algorithm dove-tails with the EV-scale, real-time algorithm. For
example, this algorithm does not determine which EVs should be used
for charging, but the real-time algorithm can make these decisions to
maximize efficiency.


Features still to incorporate:
 - Differential efficiency at different charging rates
 - Limits on charging within each site branch
