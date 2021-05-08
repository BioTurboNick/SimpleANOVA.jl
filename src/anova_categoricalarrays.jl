using .CategoricalArrays

import Base.:*

*(x::CategoricalValue, y) = unwrap(x) * y
