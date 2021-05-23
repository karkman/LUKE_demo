#!/bin/bash

esearch -db sra -query $1 |efetch --format runinfo |awk -F "," 'NR>1 {print $0}' |cut -d "," -f 1,30|sed 's/_/,/g'
