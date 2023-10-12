#ifndef INHERITANCE_H
#define INHERITANCE_H

// ͷ�ļ�����
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

// �����ض�������

// ����λ��
#define NUM 32
// ��������
#define GENES 10
// ��������
using inheri_type = std::bitset<NUM>;
// �����ֵ����
using value_type = double;
// ��ֵ��
using value_map_type = std::map<unsigned, value_type>;

#define pdTRUE 0
#define pdFALSE 1

// ������
class Inheritance
{
public:
    /* ��ʼ������ */
    Inheritance() : data(0), value(0) {}
    Inheritance(Inheritance &) = default;
    Inheritance(Inheritance &&) = default;
    ~Inheritance() = default;

    /* ��������� */
    Inheritance &operator=(Inheritance &);
    Inheritance &operator=(Inheritance &&);

    bool operator==(Inheritance &);
    bool operator==(Inheritance &&);
    bool operator!=(Inheritance &);
    bool operator!=(Inheritance &&);

    // �����ʼ�� inheri_type ����
    Inheritance(inheri_type in) : data(in), value(0)
    {
        updateValue(in);
    }

    // ���� Value
    void updateValue(inheri_type in);

    // �������»�����
    void setInheriCode(inheri_type in)
    {
        this->data = in;
        updateValue(in);
    }

    // uint8_t��uint16_t��uint32_t��uint64_t��inheri_type �Ļ������ģ��
    template <typename T>
    uint8_t setInheriCode(T in)
    {
        // ���ֱ���� inheri_type
        if (typeid(T) == typeid(inheri_type))
        {
            setInheriCode(in);
            return pdTRUE;
        }

        // ��������� uint8_t��uint16_t��uint32_t��uint64_t
        if (typeid(T) == typeid(uint32_t) || typeid(T) == typeid(unsigned) || typeid(T) == typeid(uint16_t) || typeid(T) == typeid(unsigned long) || typeid(T) == typeid(uint8_t) || typeid(T) == typeid(unsigned char) || typeid(T) == typeid(unsigned short) || typeid(T) == typeid(uint64_t))
        {
            unsigned bit_num = sizeof(T) * 8;

            // ���������Ƚ�Сʱ
            if (NUM < bit_num)
            {
                inheri_type tmp;
                for (int i = 0; i < NUM; ++i)
                {
                    tmp.set(i, static_cast<bool>(in & (1 << i)));
                }

                setInheriCode(tmp);
            }
            // ���������Ƚϴ�ʱ
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

        // ���ɶҲ����ֱ���׳�����
        std::string typeinfo(typeid(T).name());
        throw std::runtime_error("Fatal Error: Type Error (Type name: " + typeinfo + ")\n");
    }

    // ���û�����
    void reset(void)
    {
        this->data = 0;
        this->value = 0;
    }

    // ���򽻲溯��
    std::pair<inheri_type, inheri_type> cross(Inheritance &in);
    // ���򽻲溯��(ָ���)
    std::pair<inheri_type, inheri_type> cross(Inheritance *in);

    // ����ͻ�亯��
    inheri_type variate(unsigned count);

    // ����������жϺ���
    bool judge(void);

    // �������Է��غ���
    inheri_type getGene(void)
    {
        return data;
    }

    value_type getValue(void)
    {
        return value;
    }

private:
    // ��������
    inheri_type data;
    // ��ֵ����
    value_type value;
};

// ��������

// ˢ��α�����
void randomFlush(unsigned count);
// ����ɸѡ����
void train(void);
// �����ֵ���㺯��
template <std::size_t N>
value_type calculateValue(std::bitset<N> &in);

#endif