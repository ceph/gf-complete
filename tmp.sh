for i in 5 10 ; do
  sed 's/1 }/'$i' }/' tmp-time-test.sh > tmp2.sh
  sh tmp2.sh 4 LOG - - >> tmp-$i-out.txt
  sh tmp2.sh 4 TABLE - - >> tmp-$i-out.txt
  sh tmp2.sh 4 TABLE SINGLE,SSE - >> tmp-$i-out.txt
  sh tmp2.sh 8 LOG - - >> tmp-$i-out.txt
  sh tmp2.sh 8 TABLE - - >> tmp-$i-out.txt
  sh tmp2.sh 8 SPLIT 8 4 SSE - >> tmp-$i-out.txt
  sh tmp2.sh 16 LOG - - >> tmp-$i-out.txt
  sh tmp2.sh 16 SPLIT 16 4 SSE,STDMAP - >> tmp-$i-out.txt
  sh tmp2.sh 16 SPLIT 16 4 SSE,ALTMAP - >> tmp-$i-out.txt
  sh tmp2.sh 32 SPLIT 8 8 - - >> tmp-$i-out.txt
  sh tmp2.sh 32 SPLIT 32 4 SSE,STDMAP - >> tmp-$i-out.txt
  sh tmp2.sh 32 SPLIT 32 4 SSE,ALTMAP - >> tmp-$i-out.txt
done
