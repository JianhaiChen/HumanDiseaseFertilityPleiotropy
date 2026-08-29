#!/bin/bash
python "$3" --rg "$1.sumstats.gz","$2.sumstats.gz" --ref-ld-chr "$4" --w-ld-chr "$4" --out "$5" >/dev/null 2>&1
