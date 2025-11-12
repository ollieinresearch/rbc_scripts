#!/bin/bash

# Parameters given from folder organization:
CURRENT_DIR=$(pwd)
RES_GAM_DIR=$(echo "$CURRENT_DIR" | grep -oP '\d+_Gam\d+')

# Rayleigh Number
RA=$(echo "$CURRENT_DIR" | grep -oP 'ra[0-9eE\+\.-]+')

# Exponent of 10 for Pr. (ie if Pr=1, PR_EXP=0, and Pr=0.1 -> PR_EXP=-1)
PR_EXP=$(echo "$CURRENT_DIR" | grep -oP 'pr[0-9eE\+\.-]+')

# Vertical resolution
RES=$(echo "$RES_GAM_DIR" | grep -oP '^\d+')

# Aspect ratio
GAM=$(echo "$RES_GAM_DIR" | grep -oP '(?<=_Gam)\d+')