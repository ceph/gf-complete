echo "SHIFT" `gf_time 32 M 0 10240 10240 SHIFT  - - | tail -n 1`
echo "GROUP 2 4" `gf_time 32 M 0 10240 10240 GROUP 2 4  - - | tail -n 1`
echo "GROUP 3 4" `gf_time 32 M 0 10240 10240 GROUP 3 4  - - | tail -n 1`
echo "GROUP 4 4" `gf_time 32 M 0 10240 10240 GROUP 4 4  - - | tail -n 1`
echo "GROUP 2 8" `gf_time 32 M 0 10240 10240 GROUP 2 8  - - | tail -n 1`
echo "GROUP 3 8" `gf_time 32 M 0 10240 10240 GROUP 3 8  - - | tail -n 1`
echo "GROUP 4 8" `gf_time 32 M 0 10240 10240 GROUP 4 8  - - | tail -n 1`
echo "GROUP 2 2" `gf_time 32 M 0 10240 10240 GROUP 2 2  - - | tail -n 1`
echo "GROUP 3 3" `gf_time 32 M 0 10240 10240 GROUP 3 3  - - | tail -n 1`
echo "BYTWO_p" `gf_time 32 M 0 10240 10240 BYTWO_p  - - | tail -n 1`
echo "BYTWO_b" `gf_time 32 M 0 10240 10240 BYTWO_b  - - | tail -n 1`
echo "SPLIT 32 2" `gf_time 32 M 0 10240 10240 SPLIT 32 2  - - | tail -n 1`
echo "SPLIT 32 4" `gf_time 32 M 0 10240 10240 SPLIT 32 4  - - | tail -n 1`
echo "SPLIT 32 8" `gf_time 32 M 0 10240 10240 SPLIT 32 8  - - | tail -n 1`
echo "SPLIT 8 8" `gf_time 32 M 0 10240 10240 SPLIT 8 8  - - | tail -n 1`
echo "COMPOSITE 2 16 -" `gf_time 32 M 0 10240 10240 COMPOSITE 2 16 -  - - | tail -n 1`
