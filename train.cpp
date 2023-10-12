// �Ŵ��㷨
#include "train.h"

// �Ƚ����������
template <std::size_t N>
struct compBitset
{
    bool operator()(std::bitset<N> &lbs, std::bitset<N> &rbs)
    {
        // �Ƚ��㷨
        value_type value_left = calculateValue<N>(lbs);
        value_type value_right = calculateValue<N>(rbs);

        return value_left > value_right;
    }
};

// ȫ��������
namespace
{
    // ȫ�����������
    std::default_random_engine randomEngine;
    std::uniform_int_distribution<unsigned> randomVariate(0, NUM - 1);
    std::uniform_real_distribution<double> randomExample(0, 1);

    // ����̽����
    double exploreRate = 0.3;

    // ѵ������
    unsigned epoch = 10;
    const unsigned start = epoch;

    // ��ֵ����
    value_type TopLine = 1000;

    // ����λ��Ӧ��ֵ��(������32λ����)
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
    // // ��ʼ�����Ժ���
    // for (int i = 0; i < NUM; ++i)
    // {
    //     std::cout << "Value Map: " << i << " - " << valueMap[i] << std::endl;
    // }

    // ѵ������
    train();

    system("pause");
    return 0;
}

void train(void)
{
    // ʮ����ʼ��������������ֻ���10������
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

    // ʮ����ѻ��򱣴��� genes_best��������ʱ��������󱣴��� genes_all
    /* ʹ�� sort ��������򣬻��� unique ��ʵ��ȥ�� */
    std::vector<inheri_type> genes_best, genes_all;

    // ���ݽ����еĳ�ʼ����
    for (auto i : genes_groups)
    {
        genes_all.push_back(inheri_type(i));
        genes_best.push_back(inheri_type(i));
    }

    // ����ָ��洢 Inheritance ��ָ������
    // vec:([ptr:(Inheritance)][...][...][...][...][...][...][...][...][...])
    std::vector<std::unique_ptr<Inheritance>> vec;

    // ֻ��ʮ��������Բ��뽻��
    assert(genes_best.size() == GENES && "��ѻ���������ƥ�䣡");

    for (int i = 0; i < GENES; ++i)
    {
        inheri_type in = genes_best[i];
        vec.push_back(std::unique_ptr<Inheritance>(new Inheritance(in)));
    }

    // �����ѻ���
    genes_best.clear();

    // ��ʼ����ѵ��
    while (epoch--)
    {
        // Epoch Debug ��Ϣ
        std::cout << "/-----------------------------\n\
        Epoch: " << start - epoch
                  << std::endl;

        std::pair<inheri_type, inheri_type> genes;

        /* �����з��ܵĻ������ genes_all (std::set<inheri_type>) */
        // ������򣬿��Ժ��Լ����棬��Ϊ�б���Ĺ���
        for (int i = 0; i < 10; ++i)
        {
            for (int j = 0; j < 10; ++j)
            {
                // �����ж���������Ƿ����
                if (vec[i]->judge())
                {
                    genes_all.push_back(vec[i]->getGene());
                }

                if (vec[j]->judge())
                {
                    genes_all.push_back(vec[j]->getGene());
                }

                // ���л��򽻴�
                genes = vec[i]->cross(vec[j].get());

                // ����������װ�뼯�ϲ�����ͬʱȥ���ظ��Ļ���
                if (vec[i]->judge())
                {
                    genes_all.push_back(genes.first);
                }

                if (vec[j]->judge())
                {
                    genes_all.push_back(genes.second);
                }

                // ���λ������
                vec[i]->variate(randomVariate(randomEngine));
                vec[j]->variate(randomVariate(randomEngine));

                // �������Ļ���
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

        // vector ���������ȥ��
        std::sort(genes_all.begin(), genes_all.end(), compBitset<NUM>());
        auto pEnd = std::unique(genes_all.begin(), genes_all.end());
        genes_all.erase(pEnd, genes_all.end());

        // ������ѻ�����������ĸ���
        unsigned exp_num = static_cast<unsigned>(exploreRate * 10);
        unsigned bst_num = 10 - exp_num;

        auto tmp_num = genes_all.size();
        tmp_num = (tmp_num > bst_num) ? bst_num : tmp_num;

        std::reverse(genes_best.begin(), genes_best.end());

        // ������ѻ���
        for (unsigned i = 0; i < tmp_num; ++i)
        {
            genes_best.push_back(genes_all[i]);
            // genes_best.push_back(*back_it++);
        }

        for (unsigned i = 0, count = 0; i < exp_num;)
        {
            double rate = randomExample(randomEngine);
            std::size_t sz = genes_all.size();

            // �ٷֱ����ѡ��
            unsigned randomRst = static_cast<unsigned>(rate * sz);

            if (genes_all.size() == 0)
            {
                break;
            }

            // ע�� genes_all Ϊ������
            if (randomRst == genes_all.size())
            {
                randomRst--;
            }

            inheri_type rdm = genes_all[randomRst];

            // �����������
            auto tmp = calculateValue<NUM>(rdm);
            if (tmp < TopLine)
            {
                genes_best.push_back(rdm);
                ++i;
            }

            // ���û�з��ϵĻ���
            if (count++ > 20)
                break;
        }

        if (genes_best.size() < 10)
        {
            std::size_t sz = 10 - genes_all.size();

            // �����ѡ��Ļ�����10�����򲹳�ȫ�����
            for (std::size_t k = 0; k < sz; ++k)
            {
                genes_best.push_back(inheri_type(0x00000000));
            }
        }

        // �����������
        genes_all.clear();

        // vector ���������ȥ��
        std::sort(genes_best.begin(), genes_best.end(), compBitset<NUM>());
        auto pBEnd = std::unique(genes_best.begin(), genes_best.end());
        genes_best.erase(pBEnd, genes_best.end());

        // ������ѻ���ͬʱ��ӡ�����־
        std::cout << "Genes: [ " << std::endl;
        for (int i = 0; i < GENES; ++i)
        {
            auto code = genes_best[i];
            vec[i]->setInheriCode(code);

            std::cout << code << "\t" << vec[i]->getValue() << std::endl;
        }

        std::cout << "]\n"
                  << std::endl;

        // �����µĴ���ռ䣬�����ѻ��򣬲��Ҽ�����г�Ա
        // genes_best.clear();
    }
}

// ���򽻴�����ֱ���������������
std::pair<inheri_type, inheri_type> Inheritance::cross(Inheritance &in)
{
    // ��ˢ��α�������ֹ���㽻��
    randomFlush(randomVariate(randomEngine));

    unsigned var_pos = randomVariate(randomEngine);
    std::bitset<NUM> var_gene1(0x0000), var_gene2(0x0000);

    // ���λ���򽻴�(����������)
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

    // �������ȷ��
    this->data = var_gene1;
    in.data = var_gene2;

    // �������ߵĽ�������
    return std::make_pair(this->data, in.data);
}

// ���򽻴�����ֱ���������������
std::pair<inheri_type, inheri_type> Inheritance::cross(Inheritance *in)
{
    // ��ˢ��α�������ֹ���㽻��
    randomFlush(randomVariate(randomEngine));

    unsigned var_pos = randomVariate(randomEngine);
    std::bitset<NUM> var_gene1(0x0000), var_gene2(0x0000);

    // ���λ���򽻴�(����������)
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

    // �������ȷ��
    this->data = var_gene1;
    in->data = var_gene2;

    // �������ߵĽ�������
    return std::make_pair(this->data, in->data);
}

// ����ͻ�亯����ֱ���ڻ����Ͻ����޸�
inheri_type Inheritance::variate(unsigned count)
{
    // �������ͻ�� 50%
    if (count > NUM / 2)
    {
        count = NUM / 2;
    }

    for (unsigned i = 0; i < count; ++i)
    {
        // ���ͻ�� bit(�������λȡ��)
        unsigned bit = randomVariate(randomEngine);

        this->data.set(bit, (this->data.test(bit) ^ 1));
    }

    return this->data;
}

// ������ѡ����
bool Inheritance::judge(void)
{
    updateValue(this->data);
    return this->value < TopLine;
}

// ���»��� Value
void Inheritance::updateValue(inheri_type in)
{
    // ��������� Value
    this->value = calculateValue<NUM>(in);
}

// std::bitset<N> ���͵� Value ����ģ��
template <std::size_t N>
value_type calculateValue(std::bitset<N> &in)
{
    value_type value = 0;

    // ��ע��Խ������
    for (unsigned i = 0; i < N; ++i)
    {
        if (in.test(i))
            value += valueMap[i];
    }

    return value;
}

// ˢ�����������
void randomFlush(unsigned count)
{
    // ���ʹ������ˢ��α�����
    while (count-- < 0)
        randomVariate(randomEngine);
}

/* ���¾�Ϊ��������� */

// ��ֵ��ֵ�����
Inheritance &Inheritance::operator=(Inheritance &in)
{
    this->data = in.data;
    this->value = in.value;

    return *this;
}

// ��ֵ��ֵ�����
Inheritance &Inheritance::operator=(Inheritance &&in)
{
    this->data = in.data;
    this->value = in.value;

    return *this;
}

// ��ֵ��������
bool Inheritance::operator==(Inheritance &in)
{
    return this->data == in.data;
}

// ��ֵ��������
bool Inheritance::operator==(Inheritance &&in)
{
    return this->data == in.data;
}

// ��ֵ���������
bool Inheritance::operator!=(Inheritance &in)
{
    return this->data != in.data;
}

// ��ֵ���������
bool Inheritance::operator!=(Inheritance &&in)
{
    return this->data != in.data;
}
