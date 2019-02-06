Not(args...) = Not(args)
Not(args::Tuple) = args |> collect |> Not
