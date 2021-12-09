using Documenter, GigaFV
#https://juliadocs.github.io/Documenter.jl/stable/man/guide/index.html 
# Run this by julia --project=. --color=yes docs/make.jl in the root
# of this repository.
makedocs(
    modules = [GigaFV],
    sitename="GigaFV Documentation")