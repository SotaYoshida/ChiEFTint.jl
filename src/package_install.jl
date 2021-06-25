using Pkg

ts = ["Combinatorics",
      "LaTeXStrings",
      "LinearAlgebra",
      "Printf",
      "PyCall",
      "SpecialFunctions",
      "TimerOutputs",
      "ThreadPools",
      "WignerSymbols"]
for tmp in ts
    Pkg.add(tmp)
end
