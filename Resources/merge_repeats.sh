
gunzip -c $1 > $1.temp
gunzip -c $2 > $2.temp

cat $1.temp $2.temp > $3

rm $1 $2 $1.temp $2.temp
