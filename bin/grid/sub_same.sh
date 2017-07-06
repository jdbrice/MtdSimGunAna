# @Author: Daniel
# @Date:   2017-05-12 10:32:53
# @Last Modified by:   Daniel
# @Last Modified time: 2017-05-12 11:05:33


source ~/.bash_profile
wd=`pwd`
echo "pwd: " `pwd`
echo "config: $wd/$1"
cmd='Arguments='"$wd/$1"' --jobIndex=$(Process)'
condor_submit grid/same.condor.submit.sh -append "$cmd"
# echo $cmd
# echo "condor_submit grid/same.condor.submit.sh -append $cmd"
# condor_submit grid/same.condor.submit.sh -append $cmd
# for i in {0..49} 
# do
# 	echo "condor_submit grid/same.condor.submit.sh -append \"$cmd$i\""
# 	condor_submit grid/same.condor.submit.sh -append "$cmd$i"
# done
