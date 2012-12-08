if [ $# -lt 4 ]; then
  echo 'usage: sh tmp-test.sh w gf_specs (e.g. LOG - -)' >&2
  exit 1
fi

w=$1
shift
i=1024
while [ $i -le 1073741824 ]; do
  iter=`echo $i | awk '{ print (1073741824/$1)*10 }'`
  echo $i $iter $w $* `gf_time $w R -1 $i $iter $*`
  i=`echo $i | awk '{ print $1*2 }'`
done
