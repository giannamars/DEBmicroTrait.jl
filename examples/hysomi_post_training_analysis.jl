using DEBmicroTrait
using Serialization

const HYSOMI_POST_OUTPUT_PATH = joinpath(@__DIR__, "data", "outputs")

function run_hysomi_post_training_analysis_example(output_library_path = HYSOMI_POST_OUTPUT_PATH, hash_id = nothing)
    run_dirs = filter(name -> startswith(name, "train_"), readdir(output_library_path))
    isempty(run_dirs) && throw(ArgumentError("no training runs found in $output_library_path"))
    run_dir = hash_id === nothing ? joinpath(output_library_path, first(run_dirs)) : joinpath(output_library_path, "train_" * hash_id)
    return DEBmicroTrait.hysomi_post_training_analysis(output_library_path, split(basename(run_dir), "train_")[end])
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_post_training_analysis_example()
    println(result.run_class)
    println(result.evaluation_result)
end
