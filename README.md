# 2D-Simulation
2D simulation
5.5
今日成功本机测试openmp,发现单核跑两个小时的东西双核跑大概要一个小时二十分钟，也就是说开销还挺大（/(ㄒoㄒ)/~~），也就可以解释昨天12个核大概也就快了三四倍的结果了。另外，今天找到了两篇文章从题目看可以通过深度学习加速MC方法，但是仔细看发现针对的是伊辛模型和另一个凝聚态的模型，似乎对我没什么用。。今天的一个脑洞是猜测openmp的效果不好是因为main函数定义的变量比较多，在运算过程中不同线程想调用它们需要等。如果这样的话，我可以写一个函数把蒙卡包进去，在函数里初始化这些参数，这样就互不相干且每次运行完都可以销毁掉不占内存。什么时候可以试一下。