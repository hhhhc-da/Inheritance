// 遗传算法
#include "train.h"

// 比较运算符重载
template <std::size_t N>
struct compBitset
{
    bool operator()(std::bitset<N> &lbs, std::bitset<N> &rbs)
    {
        // 比较算法
        value_type value_left = calculateValue<N>(lbs);
        value_type value_right = calculateValue<N>(rbs);

        return value_left > value_right;
    }
};

// 全局量定义
namespace
{
    // 全局随机数引擎
    std::default_random_engine randomEngine;
    std::uniform_int_distribution<unsigned> randomVariate(0, NUM - 1);
    std::uniform_real_distribution<double> randomExample(0, 1);

    // 基因探索率
    double exploreRate = 0.3;

    // 训练次数
    unsigned epoch = 10;
    const unsigned start = epoch;

    // 价值容量
    value_type TopLine = 1000;

    // 基因位对应价值表(适用于32位基因)
    value_map_type valueMap{
        std::make_pair(0, 20),
        std::make_pair(1, 15),
        std::make_pair(2, 2),
        std::make_pair(3, 81),
        std::make_pair(4, 12),
        std::make_pair(5, 42),
        std::make_pair(6, 29),
        std::make_pair(7, 30),
        std::make_pair(8, 71),
        std::make_pair(9, 56),
        std::make_pair(10, 70),
        std::make_pair(11, 35),
        std::make_pair(12, 64),
        std::make_pair(13, 6),
        std::make_pair(14, 62),
        std::make_pair(15, 41),
        std::make_pair(16, 20),
        std::make_pair(17, 59),
        std::make_pair(18, 102),
        std::make_pair(19, 49),
        std::make_pair(20, 26),
        std::make_pair(21, 75),
        std::make_pair(22, 82),
        std::make_pair(23, 14),
        std::make_pair(24, 5),
        std::make_pair(25, 95),
        std::make_pair(26, 38),
        std::make_pair(27, 22),
        std::make_pair(28, 7),
        std::make_pair(29, 43),
        std::make_pair(30, 27),
        std::make_pair(31, 18)};
}

int main()
{
    // // 初始化测试函数
    // for (int i = 0; i < NUM; ++i)
    // {
    //     std::cout << "Value Map: " << i << " - " << valueMap[i] << std::endl;
    // }

    // 训练函数
    train();

    system("pause");
    return 0;
}

void train(void)
{
    // 十条初始基因链（适用于只存活10条基因）
    std::vector<inheri_type> genes_groups{
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000,
        0x00000000};

    // 十根最佳基因保存在 genes_best，所有临时基因排序后保存在 genes_all
    /* 使用 sort 库进行排序，还有 unique 库实现去重 */
    std::vector<inheri_type> genes_best, genes_all;

    // 传递进所有的初始基因
    for (auto i : genes_groups)
    {
        genes_all.push_back(inheri_type(i));
        genes_best.push_back(inheri_type(i));
    }

    // 智能指针存储 Inheritance 类指针向量
    // vec:([ptr:(Inheritance)][...][...][...][...][...][...][...][...][...])
    std::vector<std::unique_ptr<Inheritance>> vec;

    // 只有十条基因可以参与交叉
    assert(genes_best.size() == GENES && "最佳基因数量不匹配！");

    for (int i = 0; i < GENES; ++i)
    {
        inheri_type in = genes_best[i];
        vec.push_back(std::unique_ptr<Inheritance>(new Inheritance(in)));
    }

    // 清空最佳基因
    genes_best.clear();

    // 开始计数训练
    while (epoch--)
    {
        // Epoch Debug 信息
        std::cout << "/-----------------------------\n\
        Epoch: " << start - epoch
                  << std::endl;

        std::pair<inheri_type, inheri_type> genes;

        /* 将所有繁衍的基因放入 genes_all (std::set<inheri_type>) */
        // 两组基因，可以和自己交叉，因为有变异的过程
        for (int i = 0; i < 10; ++i)
        {
            for (int j = 0; j < 10; ++j)
            {
                // 首先判断自身基因是否符合
                if (vec[i]->judge())
                {
                    genes_all.push_back(vec[i]->getGene());
                }

                if (vec[j]->judge())
                {
                    genes_all.push_back(vec[j]->getGene());
                }

                // 进行基因交错
                genes = vec[i]->cross(vec[j].get());

                // 将两条基因装入集合并排序，同时去掉重复的基因
                if (vec[i]->judge())
                {
                    genes_all.push_back(genes.first);
                }

                if (vec[j]->judge())
                {
                    genes_all.push_back(genes.second);
                }

                // 随机位数变异
                vec[i]->variate(randomVariate(randomEngine));
                vec[j]->variate(randomVariate(randomEngine));

                // 检查变异后的基因
                if (vec[i]->judge())
                {
                    genes_all.push_back(vec[i]->getGene());
                }

                if (vec[j]->judge())
                {
                    genes_all.push_back(vec[j]->getGene());
                }
            }
        }

        // vector 进行排序和去重
        std::sort(genes_all.begin(), genes_all.end(), compBitset<NUM>());
        auto pEnd = std::unique(genes_all.begin(), genes_all.end());
        genes_all.erase(pEnd, genes_all.end());

        // 计算最佳基因和随机基因的个数
        unsigned exp_num = static_cast<unsigned>(exploreRate * 10);
        unsigned bst_num = 10 - exp_num;

        auto tmp_num = genes_all.size();
        tmp_num = (tmp_num > bst_num) ? bst_num : tmp_num;

        std::reverse(genes_best.begin(), genes_best.end());

        // 加入最佳基因
        for (unsigned i = 0; i < tmp_num; ++i)
        {
            genes_best.push_back(genes_all[i]);
            // genes_best.push_back(*back_it++);
        }

        for (unsigned i = 0, count = 0; i < exp_num;)
        {
            double rate = randomExample(randomEngine);
            std::size_t sz = genes_all.size();

            // 百分比随机选择
            unsigned randomRst = static_cast<unsigned>(rate * sz);

            if (genes_all.size() == 0)
            {
                break;
            }

            // 注意 genes_all 为零的情况
            if (randomRst == genes_all.size())
            {
                randomRst--;
            }

            inheri_type rdm = genes_all[randomRst];

            // 插入随机基因
            auto tmp = calculateValue<NUM>(rdm);
            if (tmp < TopLine)
            {
                genes_best.push_back(rdm);
                ++i;
            }

            // 如果没有符合的基因
            if (count++ > 20)
                break;
        }

        if (genes_best.size() < 10)
        {
            std::size_t sz = 10 - genes_all.size();

            // 如果遴选后的基因不足10个，则补充全零基因
            for (std::size_t k = 0; k < sz; ++k)
            {
                genes_best.push_back(inheri_type(0x00000000));
            }
        }

        // 清除所有容器
        genes_all.clear();

        // vector 进行排序和去重
        std::sort(genes_best.begin(), genes_best.end(), compBitset<NUM>());
        auto pBEnd = std::unique(genes_best.begin(), genes_best.end());
        genes_best.erase(pBEnd, genes_best.end());

        // 更新最佳基因，同时打印输出日志
        std::cout << "Genes: [ " << std::endl;
        for (int i = 0; i < GENES; ++i)
        {
            auto code = genes_best[i];
            vec[i]->setInheriCode(code);

            std::cout << code << "\t" << vec[i]->getValue() << std::endl;
        }

        std::cout << "]\n"
                  << std::endl;

        // 创建新的代码空间，清空最佳基因，并且检查所有成员
        // genes_best.clear();
    }
}

// 基因交错函数会直接作用于两侧基因
std::pair<inheri_type, inheri_type> Inheritance::cross(Inheritance &in)
{
    // 先刷新伪随机数防止定点交换
    randomFlush(randomVariate(randomEngine));

    unsigned var_pos = randomVariate(randomEngine);
    std::bitset<NUM> var_gene1(0x0000), var_gene2(0x0000);

    // 随机位基因交错(基因拆分重组)
    for (unsigned i = 0; i < var_pos; ++i)
    {
        for (unsigned i = 0; i < var_pos; ++i)
        {
            var_gene1.set(i, this->data.test(i));
            var_gene2.set(i, in.data.test(i));
        }

        for (unsigned i = var_pos; i < NUM; ++i)
        {
            var_gene1.set(i, in.data.test(i));
            var_gene2.set(i, this->data.test(i));
        }
    }

    // 基因更改确认
    this->data = var_gene1;
    in.data = var_gene2;

    // 返回两者的交错后基因
    return std::make_pair(this->data, in.data);
}

// 基因交错函数会直接作用于两侧基因
std::pair<inheri_type, inheri_type> Inheritance::cross(Inheritance *in)
{
    // 先刷新伪随机数防止定点交换
    randomFlush(randomVariate(randomEngine));

    unsigned var_pos = randomVariate(randomEngine);
    std::bitset<NUM> var_gene1(0x0000), var_gene2(0x0000);

    // 随机位基因交错(基因拆分重组)
    for (unsigned i = 0; i < var_pos; ++i)
    {
        for (unsigned i = 0; i < var_pos; ++i)
        {
            var_gene1.set(i, this->data.test(i));
            var_gene2.set(i, in->data.test(i));
        }

        for (unsigned i = var_pos; i < NUM; ++i)
        {
            var_gene1.set(i, in->data.test(i));
            var_gene2.set(i, this->data.test(i));
        }
    }

    // 基因更改确认
    this->data = var_gene1;
    in->data = var_gene2;

    // 返回两者的交错后基因
    return std::make_pair(this->data, in->data);
}

// 基因突变函数会直接在基因上进行修改
inheri_type Inheritance::variate(unsigned count)
{
    // 最多允许突变 50%
    if (count > NUM / 2)
    {
        count = NUM / 2;
    }

    for (unsigned i = 0; i < count; ++i)
    {
        // 随机突变 bit(基因变异位取反)
        unsigned bit = randomVariate(randomEngine);

        this->data.set(bit, (this->data.test(bit) ^ 1));
    }

    return this->data;
}

// 基因遴选函数
bool Inheritance::judge(void)
{
    updateValue(this->data);
    return this->value < TopLine;
}

// 更新基因 Value
void Inheritance::updateValue(inheri_type in)
{
    // 更新自身的 Value
    this->value = calculateValue<NUM>(in);
}

// std::bitset<N> 类型的 Value 计算模板
template <std::size_t N>
value_type calculateValue(std::bitset<N> &in)
{
    value_type value = 0;

    // 请注意越界问题
    for (unsigned i = 0; i < N; ++i)
    {
        if (in.test(i))
            value += valueMap[i];
    }

    return value;
}

// 刷新随机数函数
void randomFlush(unsigned count)
{
    // 多次使用引擎刷新伪随机数
    while (count-- < 0)
        randomVariate(randomEngine);
}

/* 以下均为运算符定义 */

// 左值赋值运算符
Inheritance &Inheritance::operator=(Inheritance &in)
{
    this->data = in.data;
    this->value = in.value;

    return *this;
}

// 右值赋值运算符
Inheritance &Inheritance::operator=(Inheritance &&in)
{
    this->data = in.data;
    this->value = in.value;

    return *this;
}

// 左值相等运算符
bool Inheritance::operator==(Inheritance &in)
{
    return this->data == in.data;
}

// 右值相等运算符
bool Inheritance::operator==(Inheritance &&in)
{
    return this->data == in.data;
}

// 左值不等运算符
bool Inheritance::operator!=(Inheritance &in)
{
    return this->data != in.data;
}

// 右值不等运算符
bool Inheritance::operator!=(Inheritance &&in)
{
    return this->data != in.data;
}
