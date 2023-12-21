for x in $(ls /data/reads/BIOS2936); do echo ${x%???}; zcat /data/reads/BIOS2936/$x | awk '{if (NR%4==2)  print substr($1,1,50); else if (NR%4==0) print substr($1,1,50); else print}' > ${x%???}; gzip ${x%???}; done


${x%???} trims .gz