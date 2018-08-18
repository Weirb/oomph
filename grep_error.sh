#!/bin/bash
grep "Norm of error" output.* | awk 'BEGIN{FS=":";ORS="\n"} {print $NF}'
