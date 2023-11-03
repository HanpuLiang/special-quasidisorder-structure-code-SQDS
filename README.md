# specical-quasirandom-structure-SQS

用来存放一些SQS和SQDS的脚本。

## SQS

对于脚本sqs2poscar-HPLiang.py，这个脚本没有原子数量限制。GitHub上比较经典的那个sqs2poscar有原子数量400的限制，所以搞超大胞的时候会报错。

执行方法是
`python sqs2poscar-HPLiang.py bestsqs.out`

## SQDS

这个是代替ATAT的SQS方法，需要指定关联函数。需要的文件有：超胞的POSCAR文件、团簇文件clusters.out、结构输入文件rndstr.in

需要在python文件中修改一下超胞文件名、clusters.out中各团簇对应的关联函数、各团簇的权重、可交换原子的序号的集合。

执行方法
`python SQDS-HPLiang.py`

运行成功的话，会生成一个save-best-data目录，里面会存放所有搜索到的结构，等到跑一段时间以后，可以取score最小的那个结构，也就是拟合最好的。