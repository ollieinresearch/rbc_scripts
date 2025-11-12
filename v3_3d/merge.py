from dedalus.tools.post import visit_writes, merge_data, merge_setup

# Suppose merged files are at paths in merged_paths = ["merged_set1.h5", "merged_set2.h5", ...]
final = "snapshots.h5"

# First, set up the final file with metadata from one of them
merge_setup(final, [], virtual=False)

# Then loop through writes and merge them
def merge_fn(set_path, start, count, final_file=None):
    # here set_path corresponds to a merged set file, and start/count to write ranges
    # for each processor chunk (though merged sets should be atomic per write)
    merge_data(final_file, set_path)  # or custom logic to append

visit_writes(merged_paths, merge_fn, final_file=final)
