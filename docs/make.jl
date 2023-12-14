using Documenter
using Microgrids
using DocStringExtensions

makedocs(
    sitename = "Microgrids",
    format = Documenter.HTML(sidebar_sitename=false),
    modules = [Microgrids],
    pages = ["Microgrids"    => "index.md"],
   
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
