# "boss" class should always contain a "design_points" (the dataframe of selected design points) slot,
# a "gp_params" slot (a list of the final gp parameters including signal_var, length_scale, quadratic mean etc),
# a "surrogate" slot (a function that takes the conditioning parameter and output a single value),
# a "essential" support (the boundary of the essential support).

# Construct method such as "print", "plot", "summary" for the class "boss"
