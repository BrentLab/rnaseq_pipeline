#!/usr/bin/env bash
srun --mem=5000 --cpus-per-task=1 -J interactive -p interactive --pty /bin/bash -l