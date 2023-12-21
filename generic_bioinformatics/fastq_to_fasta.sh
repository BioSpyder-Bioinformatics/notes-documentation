#!/bin/bash

# number=0
# for line in $(zcat $1);
# do
# 	number=$(expr $number + 1)
# 	remainder=$(($number % 4))
# 	if [ $remainder == 2 ]; then
# 		echo $line
# 	fi
# done


# number=0

# while IFS= read -r line;
# do
# 	number=$(expr $number + 1)
# 	remainder=$(($number % 4))
# 	if [ $remainder == 2 ]; then
# 		echo $line
# 	fi
# done < $1



number=0

cat $1 | while IFS= read -r line;
do
	number=$(expr $number + 1)
	remainder=$(($number % 4))
	if [ $remainder == 2 ]; then
		echo $line
	fi
done


