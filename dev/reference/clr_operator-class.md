# Implicit centred clr matrix

The centred clr matrix is `B - u 1' - 1 c'` with `B` sparse, so `irlba`
only needs its products with vectors. `scale` holds the divisor of each
column.
