#  2026/5/16

今天修改了particle.h和particle.cpp文件，主要是将二维改到了3D，之后再改别的
3D\_Hard\_Rod/

│

├── .git/

├── .gitignore

├── README.md

│

├── src/               # 所有.cpp

│   ├── main.cpp

│   ├── system.cpp

│   ├── particles.cpp

│   └── random.cpp

│

├── include/           # 所有.h

│   ├── system.h

│   ├── particles.h

│   └── random.h

│

├── input/             # 输入文件

│   └── initial.txt

│

├── output/            # 计算结果（Git忽略）

│

├── docs/              # 算法说明、论文、推导

│

├── 3D.sln

├── 3D.vcxproj

└── 3D.vcxproj.filters

# 2026/7/7

今天要做的

先完成 `debug.dump` 的输出闭环（今天）

​	今天知道了为什么输出的文件那么大，肯定是输出文件写的有问题，具体问题是类型三的粒子嵌套在了类型二里面

​	！！**得到的经验：以后出现巨大的文件一定要来检查一下输出**

​	新增 `OutputVTK()`，用 VTK 的 `LINES` 或网格来表示双曲面。这样既能减少数据量，又比点云更直观，也更适合作为论文中的可视化结果。

安装ParaView 查看VTK代码，

再整理代码结构（今天）

最后增加科研分析功能（明天）

# 2026/7/8

使用paraview

File-> Open ->scene.pvd
