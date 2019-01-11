#!/bin/sh

count=360
while [ $count -ge 1 ] ; do
   echo $count
   
   if   [ $count -ge 100 ]; then
      mkdir $count
   elif [ $count -ge 10 ]; then
      mkdir 0$count
   else
      mkdir 00$count
   fi
   
   count=`expr $count - 1`
done
