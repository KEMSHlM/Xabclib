# Memo
> [Xabclib](http://www.abc-lib.org/Xabclib/Release/Readme-v1.03.pdf)

## 1. Package の作り方を学ぶ．
CMakeでライブラリを作成する際に，foo-config.cmakeを作成する必要がある．
foo-config.cmakeを手で作成することもできるが，install(EXPORT)コマンドを利用して，自動生成することもできる．

### メリット
- FindFoo.cmakeを作成する必要がない．
- クライアントのCMakeLists.txt でfind_package(foo)した後に，target_link_libraries(foo)と書くだけで，自作のライブラリをリンクできる．

### 手順
