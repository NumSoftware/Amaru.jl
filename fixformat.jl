#!/bin/bash
find -iname '*.jl' -execdir sed -i 's/\t/    /g' {} \;
find -iname '*.jl' -execdir sed -i 's/[ \t]*$//g' {} \;
