#!/bin/bash
while getopts "N:e:c" arg #选项后面的冒号表示该选项需要参数，Nbc表示用形如-b要传入的参数
do
    case $arg in #对参数列表进行一个遍历处理
        N)
            echo "N's arg:$OPTARG" #参数存在$OPTARG中
            N_val=$OPTARG
            ;;
        e)
            echo "e's arg:$OPTARG"
            e_val=$OPTARG
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

./eigfun --add-param '{"Mesh":{"N1":'$N_val',"N2":'$N_val',"N3":'$N_val'},"Eigen":{"EigenFunIdx":  '$e_val'}}'

