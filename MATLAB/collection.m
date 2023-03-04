clc
clear
close all
[filename, filepath]=uigetfile('*.*','请选择文件');

dir_colections = [filepath 'collections'];%定义我要所收集的函数存放的文件夹。pwd表示当前目录
if exist(dir_colections,'dir') == 7%判断该文件夹是否存在，存在删除及其内容
rmdir(dir_colections,'s');   
end
mkdir(dir_colections);%建立用于存放文件的文件夹
[fList,pList] = matlab.codetools.requiredFilesAndProducts(filename);%寻找文件的依赖文件
num_files = length(fList);%依赖文件数目

for k = 1:num_files
    file_path = char(fList(k));%寻找每个依赖文件路径并字符化
    [pathstr,f_name,ext]=fileparts(file_path);%依赖文件名
    copyfile(file_path,[dir_colections '\' f_name ext])%将依赖文件拷贝到设定的文件夹下
end