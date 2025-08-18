#!/bin/bash

# List of directories to process (relative to BBUFFER)
DIRS=(
  "ollie_rb_data/3d/ra1e6/pr3e-3/144_Gam2_2"
  "ollie_rb_data/3d/ra1e6/pr1e-3/200_Gam2_2"
  "ollie_rb_data/3d/ra1e6/pr3e-4/200_Gam2_2"
)

# Subdirectories to copy
SUBDIRS=(
  "analysis"
  "field_analysis"
  "outputs"
)

for DIR in "${DIRS[@]}"; do
  BBDIR="$BBUFFER/$DIR"
  HDIR="$HOME/$DIR"

  mkdir -p "$HDIR/state"

  RECENT=$(find "$BBDIR/state/." -maxdepth 1 -type f -exec basename {} \; | sort -V | tail -n 1)

  cp -r "$BBDIR/state/$RECENT" "$HDIR/state"

  for SUBDIR in "${SUBDIRS[@]}"; do
    cp -r "$BBDIR/$SUBDIR" "$HDIR"
  done
done
