#include <iostream>
#include <limits>
int main() {
  using localTestType = long double;
  // 二进制尾数的位数（包括隐式位）
  std::cout << "二进制尾数位数: " << std::numeric_limits<localTestType>::digits
            << std::endl;
  // 十进制下确保所有值唯一表示的位数（即有效数字位数）
  std::cout << "十进制有效位数（digits10）: "
            << std::numeric_limits<localTestType>::digits10 << std::endl;
  // 转换为字符串再转回时不丢失精度的最大十进制位数
  std::cout << "最大十进制位数（max_digits10）: "
            << std::numeric_limits<localTestType>::max_digits10 << std::endl;
}