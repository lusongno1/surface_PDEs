#!/bin/bash
while getopts "N:E:c" arg #选项后面的冒号表示该选项需要参数，Nbc表示用形如-b要传入的参数
do
    case $arg in #对参数列表进行一个遍历处理
        N)
         #   echo "a's arg:$OPTARG" #参数存在$OPTARG中
            N_val=$OPTARG
            ;;
        E)
          #  echo "b's arg:$OPTARG"
	    exe=$OPTARG

            ;;
        c)
            echo "c's arg:$OPTARG"
            ;;
        ?)  #当有不认识的选项的时候arg为?
            echo "unkonw argument"
            exit 1
        ;;
    esac
done

./$exe --add-param '{"Mesh":{"N1":'$N_val',"N2":'$N_val',"N3":'$N_val'}}'

