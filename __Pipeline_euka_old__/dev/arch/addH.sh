awk  '/^>/ {gsub(/.aa(sta)?$/,"",'$2');printf($0 "_" '$2' "\n");next;} {print}' $1 > $1_h
