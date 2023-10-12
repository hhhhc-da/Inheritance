#ifndef INHERITANCE_H
#define INHERITANCE_H

// 头文件声明
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <utility>
#include <cstddef>
#include <random>
#include <bitset>
#include <exception>
#include <typeinfo>
#include <queue>
#include <memory>
#include <assert.h>
#include <algorithm>

// 类型重定义声明

// 基因位数
#define NUM 32
// 基因存活数
#define GENES 10
// 基因类型
using inheri_type = std::bitset<NUM>;
// 基因价值类型
using value_type = double;
// 价值表
using value_map_type = std::map<unsigned, value_type>;

#define pdTRUE 0
#define pdFALSE 1

// 类声明
class Inheritance
{
public:
    /* 初始化函数 */
    Inheritance() : data(0), value(0) {}
    Inheritance(Inheritance &) = default;
    Inheritance(Inheritance &&) = default;
    ~Inheritance() = default;

    /* 运算符重载 */
    Inheritance &operator=(Inheritance &);
    Inheritance &operator=(Inheritance &&);

    bool operator==(Inheritance &);
    bool operator==(Inheritance &&);
    bool operator!=(Inheritance &);
    bool operator!=(Inheritance &&);

    // 特殊初始化 inheri_type 类型
    Inheritance(inheri_type in) : data(in), value(0)
    {
        updateValue(in);
    }

    // 更新 Value
    void updateValue(inheri_type in);

    // 基础更新基因函数
    void setInheriCode(inheri_type in)
    {
        this->data = in;
        updateValue(in);
    }

    // uint8_t、uint16_t、uint32_t、uint64_t、inheri_type 的基因更新模板
    template <typename T>
    uint8_t setInheriCode(T in)
    {
        // 如果直接是 inheri_type
        if (typeid(T) == typeid(inheri_type))
        {
            setInheriCode(in);
            return pdTRUE;
        }

        // 如果类型是 uint8_t、uint16_t、uint32_t、uint64_t
        if (typeid(T) == typeid(uint32_t) || typeid(T) == typeid(unsigned) || typeid(T) == typeid(uint16_t) || typeid(T) == typeid(unsigned long) || typeid(T) == typeid(uint8_t) || typeid(T) == typeid(unsigned char) || typeid(T) == typeid(unsigned short) || typeid(T) == typeid(uint64_t))
        {
            unsigned bit_num = sizeof(T) * 8;

            // 基因容量比较小时
            if (NUM < bit_num)
            {
                inheri_type tmp;
                for (int i = 0; i < NUM; ++i)
                {
                    tmp.set(i, static_cast<bool>(in & (1 << i)));
                }

                setInheriCode(tmp);
            }
            // 基因容量比较大时
            else
            {
                inheri_type tmp;
                for (int i = 0; i < bit_num; ++i)
                {
                    tmp.set(i, static_cast<bool>(in & (1 << i)));
                }

                setInheriCode(tmp);
            }
            return pdTRUE;
        }

        // 如果啥也不是直接抛出错误
        std::string typeinfo(typeid(T).name());
        throw std::runtime_error("Fatal Error: Type Error (Type name: " + typeinfo + ")\n");
    }

    // 重置基因函数
    void reset(void)
    {
        this->data = 0;
        this->value = 0;
    }

    // 基因交叉函数
    std::pair<inheri_type, inheri_type> cross(Inheritance &in);
    // 基因交叉函数(指针版)
    std::pair<inheri_type, inheri_type> cross(Inheritance *in);

    // 基因突变函数
    inheri_type variate(unsigned count);

    // 基因合理性判断函数
    bool judge(void);

    // 基因属性返回函数
    inheri_type getGene(void)
    {
        return data;
    }

    value_type getValue(void)
    {
        return value;
    }

private:
    // 基因数据
    inheri_type data;
    // 价值总数
    value_type value;
};

// 函数声明

// 刷新伪随机数
void randomFlush(unsigned count);
// 繁衍筛选函数
void train(void);
// 基因价值计算函数
template <std::size_t N>
value_type calculateValue(std::bitset<N> &in);

#endif