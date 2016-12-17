# gmp-sm2
C library for ECC-SM2 public encrypt

### 编译说明
1. 下载并编译对应平台的[libgmp](https://gmplib.org/)源码
2. 将编译好的静态库放入 lib目录下
3. 进入Debug 执行make 编译
有两种 一种是编译出静态库可以被调用
另外是编译出test 在平台下进行加密测试
