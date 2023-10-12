cd /Users/aniket/Documents/MarlinSim/03_code/therml/3d
julia --project=./therml_environment
using PackageCompiler
create_app("therml_environment","therML")