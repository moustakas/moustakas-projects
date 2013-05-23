#!/bin/tcsh
echo "redmapper_isedfit, first=0, last=0, /ised, /cl" | /usr/bin/nohup idl > & ~/archive/projects/redmapper/logs/chunk0.log &
echo "redmapper_isedfit, first=1, last=1, /ised, /cl" | /usr/bin/nohup idl > & ~/archive/projects/redmapper/logs/chunk1.log &
echo "redmapper_isedfit, first=2, last=2, /ised, /cl" | /usr/bin/nohup idl > & ~/archive/projects/redmapper/logs/chunk2.log &

echo "redmapper_isedfit, first=3, last=5, /ised, /cl" | /usr/bin/nohup idl > & ~/archive/projects/redmapper/logs/chunk3-5.log &
echo "redmapper_isedfit, first=6, last=7, /ised, /cl" | /usr/bin/nohup idl > & ~/archive/projects/redmapper/logs/chunk6-7.log &
