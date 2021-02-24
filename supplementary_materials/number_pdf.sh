#!/bin/bash

# requires enscript and pdftk-java, both installed through brew
input="$1"
output="${1%.pdf}-header.pdf"
pagenum=$(pdftk "$input" dump_data | grep "NumberOfPages" | cut -d":" -f2)
enscript -L1 --header='||Page $% of $=' --output - < <(for i in $(seq "$pagenum"); do echo; done) | ps2pdf - | pdftk "$input" multistamp - output $output
