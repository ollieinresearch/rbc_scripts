# Storage, backups and quotas

## Three places

| Space | Path | Backed up | Use |
|---|---|---|---|
| {term}`scratch` | `$SCRATCH` (e.g. `/scratch/ollie`) | no | where runs execute: fastest, largest quota |
| {term}`home` | `$HOME` | yes | the second copy of every run's results |
| {term}`project space` | `$PROJECT`, or `/project/def-goluskin/<user>` on Fir | yes | the group's copy, shared with David |

Scratch has no backups and may be purged of old files; see each cluster's page and
[Storage and file management](https://docs.alliancecan.ca/wiki/Storage_and_file_management).
Anything that matters must also be in home or project space.

## What is backed up

After every post-processing job, `backup_and_prune` mirrors every run under `DATA_ROOT`
into each of `BACKUP_DIRS`, keeping the path below `$SCRATCH`:

```text
/scratch/ollie/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2
  -> $HOME/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2
  -> /project/def-goluskin/ollie/ollie_rb_data/3d/ra1e5/pr1e0/48_Gam2_2
```

Copied: `analysis/`, `horizontal_analysis/`, `spectra/`, `outputs/`, `params.conf`, the
job script, the job logs, `stop_reason.txt`, and the `BACKUP_STATES` (2) newest readable
checkpoints. Not copied: snapshots, renders and older checkpoints.

Then:

1. older checkpoints are removed from each backup, so each backup holds two;
2. only if every copy succeeded, checkpoints beyond the `KEEP_STATES` (4) newest are
   deleted from scratch.

`collect_results` also copies `data/results.csv` (and `results.parquet`) into each backup,
next to the run tree (e.g. `$HOME/ollie_rb_data/3d/results.csv`).

Copies use `rsync`, so only new or changed files move; the first backup of a large tree
takes a while, later ones are quick.

## Quotas

Every space has limits on both bytes and number of files. Check often:

```bash
diskusage_report
```

Every job log starts with its output. A job that fills a quota dies in the middle of a
write with an HDF5 error that does not mention space (for example about a byte count
overflowing); see {doc}`../troubleshooting`. What uses the most:

| Output | Size | Kept |
|---|---|---|
| checkpoints | a few GB each at production resolution | 4 on scratch, 2 per backup |
| snapshots | 30 writes of 3 fields per set: by far the largest | until you delete them |
| process-file folders (virtual mode) | one file per rank per set: counts against the file quota | as above |
| analysis, spectra, outputs | small | always |

Delete snapshots once their movies and spectra are made; nothing deletes them
automatically.

## `results.csv` stays private

`data/results.csv` and `data/results.parquet` are in `.gitignore`: until the paper is
published they are kept only on the clusters and in the backups, not in the public
repository. Remove those lines from `.gitignore` to publish them.
