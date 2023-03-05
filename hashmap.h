//
// Created by Alexander Nekrasov on 20.01.2023.
//

// https://web.archive.org/web/20160304094014/http://www.arl.wustl.edu/~sailesh/download_files/Limited_Edition/hash/Dynamic%20Perfect%20Hashing-%20Upper%20and%20Lower%20Bounds.pdf

#pragma once
#ifndef LOCAL
//#pragma GCC optimize("O3,unroll-loops,no-stack-protector,inline,fast-math")
//#pragma GCC target("tune=native")
//#pragma GCC target("avx,avx2")
#endif

#include <cassert>
#include <functional>
#include <forward_list>
#include <random>
#include <chrono>
#include <list>
#include <stdexcept>
#include <memory>

namespace {
    const double HASH_MAP_C = 0.8;
    const size_t MAX_OFFSET = 3;

    const size_t PRIMES[29] = {17ul, 53ul, 97ul, 193ul, 389ul, 769ul, 1543ul, 3079ul, 6151ul, 12289ul, 24593ul, 49157ul,
                               98317ul, 196613ul, 393241ul, 786433ul, 1572869ul, 3145739ul, 6291469ul, 12582917ul,
                               25165843ul, 50331653ul, 100663319ul, 201326611ul, 402653189ul, 805306457ul, 1610612741ul,
                               3221225473ul, 4294967291ul};
//    const size_t PRIMES[15] = {17ul, 97ul, 389ul, 1543ul, 6151ul, 24593ul, 98317ul, 393241ul, 1572869ul, 6291469ul,
//                                25165843ul, 100663319ul, 402653189ul, 1610612741ul, 4294967291ul};

    size_t get_nearest_prime(size_t x) {
        int pos = 0;
        while (PRIMES[pos] < x) {
            ++pos;
        }
        return PRIMES[pos];
//        return *std::lower_bound(std::begin(PRIMES), std::end(PRIMES), x);
    }
}

template<class KeyType, class ValueType, class Hash = std::hash<KeyType>>
class HashMap {
private:
    struct Position {
        typename std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>::iterator iter_;
    };
public:
    explicit HashMap(Hash hasher = Hash());

    template<typename Iterator>
    HashMap(Iterator first, Iterator last, Hash hasher = Hash());

    HashMap(const std::initializer_list<std::pair<const KeyType, ValueType>> &list, Hash hasher = Hash());

    HashMap(const HashMap &other);

    HashMap &operator=(const HashMap &other);

    bool empty() const;

    size_t size() const;

    const Hash &hash_function() const;

    void insert(const std::pair<KeyType, ValueType> &key_value);

    bool contains(const KeyType &key);

    class const_iterator { // NOLINT
    public:
        explicit const_iterator(typename std::list<Position>::const_iterator iter);

        const_iterator();

        bool operator!=(const const_iterator &other) const;

        bool operator==(const const_iterator &other) const;

        const_iterator &operator++();

        const_iterator operator++(int);

        const std::pair<const KeyType, ValueType> &operator*() const;

        const std::pair<const KeyType, ValueType> *operator->() const;

    private:
        typename std::list<Position>::const_iterator iter_;
    };

    class iterator { // NOLINT
    public:
        explicit iterator(typename std::list<Position>::iterator iter);

        iterator();

        bool operator!=(const iterator &other) const;

        bool operator==(const iterator &other) const;

        iterator &operator++();

        iterator operator++(int);


        std::pair<const KeyType, ValueType> &operator*() const;

        std::pair<const KeyType, ValueType> *operator->() const;

    private:
        typename std::list<Position>::iterator iter_;
    };

    iterator begin();

    const_iterator begin() const;

    iterator end();

    const_iterator end() const;

    iterator find(const KeyType &key);

    const_iterator find(const KeyType &key) const;

    const ValueType &at(const KeyType &key) const;

    ValueType &operator[](const KeyType &key);

//    const ValueType &operator[](const KeyType &key) const;

    void erase(const KeyType &key);

    void swap(HashMap<KeyType, ValueType, Hash> &other);

    void clear();

private:

    static size_t S(size_t m); // NOLINT

    size_t max_sum_s() const;

    void rehash_all(const std::vector<std::pair<const KeyType, ValueType>> &inserted);

    std::unique_ptr<const_iterator> lookup_const(const KeyType &key) const;

    std::unique_ptr<iterator> lookup(const KeyType &key);


    void insert(const KeyType &key, const ValueType &value);

    struct SingleHashDataCell {
        std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>> values;
    };

    struct HashFunction {
        size_t y, k, b;

        void randomize();

        size_t operator()(size_t x, size_t p) const;
    };

    void
    fix_j(std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>> &l_j,
          std::vector<size_t> &hashes, size_t j);


    size_t m_;
    Hash hasher_;
    size_t count_;
    size_t size_;
    HashFunction hash_global_;
    std::vector<HashFunction> hash_local_;
    std::vector<std::vector<SingleHashDataCell>> t_;
    std::vector<size_t> b_;
    std::vector<size_t> m_local_;
    std::vector<size_t> s_local_;
    size_t sum_s_;
    std::list<Position> positions_;

    static std::mt19937_64 rng; // NOLINT
};

//template<typename KeyType, typename ValueType, typename Hash> std::mt19937_64 HashMap<KeyType, ValueType, Hash>::rng(
//        std::chrono::high_resolution_clock::now().time_since_epoch().count());
template<typename KeyType, typename ValueType, typename Hash> std::mt19937_64 HashMap<KeyType, ValueType, Hash>::rng(
        179);

/**********************************************************************
 *              Constructors and assignment operators                 *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::HashMap(Hash hasher) : HashMap({}, hasher) {
}

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::HashMap(const HashMap<KeyType, ValueType, Hash> &other)
        : HashMap(other.begin(), other.end(), other.hasher_) {
}

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::HashMap(const std::initializer_list<std::pair<const KeyType, ValueType>> &list,
                                           Hash hasher)
        : HashMap(list.begin(), list.end(), hasher) {
}

template<class KeyType, class ValueType, class Hash>
template<typename Iterator>
HashMap<KeyType, ValueType, Hash>::HashMap(Iterator first, Iterator last, Hash hasher) : hasher_(hasher) {
    hash_global_ = {1};
    size_ = 0;
    m_ = 0;
    while (first != last) {
        insert(first->first, first->second);
        ++first;
    }
}

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash> &
HashMap<KeyType, ValueType, Hash>::operator=(const HashMap<KeyType, ValueType, Hash> &other) {
    HashMap<KeyType, ValueType, Hash> tmp(other);
    swap(tmp);
    return *this;
}

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::swap(HashMap<KeyType, ValueType, Hash> &other) {
    std::swap(m_, other.m_);
    std::swap(hasher_, other.hasher_);
    std::swap(count_, other.count_);
    std::swap(size_, other.size_);
    std::swap(hash_global_, other.hash_global_);
    std::swap(hash_local_, other.hash_local_);
    std::swap(t_, other.t_);
    std::swap(b_, other.b_);
    std::swap(m_local_, other.m_local_);
    std::swap(s_local_, other.s_local_);
    std::swap(sum_s_, other.sum_s_);
    std::swap(positions_, other.positions_);
}

/**********************************************************************
 *                         const_iterator                             *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::const_iterator::const_iterator()
        : iter_() {
}


template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::const_iterator::const_iterator(typename std::list<Position>::const_iterator iter)
        : iter_(iter) {
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::const_iterator
HashMap<KeyType, ValueType, Hash>::find(const KeyType &key) const {
    auto it = lookup_const(key);
    if (it == nullptr) {
        return end();
    }
    return *it;
}

template<class KeyType, class ValueType, class Hash>
const std::pair<const KeyType, ValueType> *HashMap<KeyType, ValueType, Hash>::const_iterator::operator->() const {
    return &iter_->iter_->first;
}

template<class KeyType, class ValueType, class Hash>
const std::pair<const KeyType, ValueType> &HashMap<KeyType, ValueType, Hash>::const_iterator::operator*() const {
    return iter_->iter_->first;
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::const_iterator
HashMap<KeyType, ValueType, Hash>::const_iterator::operator++(int) {
    auto ret = *this;
    operator++();
    return ret;
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::const_iterator &
HashMap<KeyType, ValueType, Hash>::const_iterator::operator++() {
    ++iter_;
    return *this;
}

template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::const_iterator::operator==(const HashMap::const_iterator &other) const {
    return iter_ == other.iter_;
}

template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::const_iterator::operator!=(const HashMap::const_iterator &other) const {
    return !(*this == other); // NOLINT
}

/**********************************************************************
 *                             iterator                               *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::iterator::iterator() : iter_() {
}

template<class KeyType, class ValueType, class Hash>
HashMap<KeyType, ValueType, Hash>::iterator::iterator(typename std::list<Position>::iterator iter)
        : iter_(iter) {
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::iterator HashMap<KeyType, ValueType, Hash>::find(const KeyType &key) {
    auto it = lookup(key);
    if (it == nullptr) {
        return end();
    }
    return *it;
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::iterator HashMap<KeyType, ValueType, Hash>::iterator::operator++(int) {
    auto ret = *this;
    operator++();
    return ret;
}


template<class KeyType, class ValueType, class Hash>
std::pair<const KeyType, ValueType> &HashMap<KeyType, ValueType, Hash>::iterator::operator*() const {
    return iter_->iter_->first;
}

template<class KeyType, class ValueType, class Hash>
std::pair<const KeyType, ValueType> *HashMap<KeyType, ValueType, Hash>::iterator::operator->() const {
    return &(operator*());
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::iterator &HashMap<KeyType, ValueType, Hash>::iterator::operator++() {
    ++iter_;
    return *this;
}

template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::iterator::operator==(const HashMap::iterator &other) const {
    return iter_ == other.iter_;
}

template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::iterator::operator!=(const HashMap::iterator &other) const {
    return !(*this == other); // NOLINT
}

/**********************************************************************
 *                       access functions                             *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
ValueType &HashMap<KeyType, ValueType, Hash>::operator[](const KeyType &key) {
    auto it = lookup(key);
    if (it == nullptr) {
        ValueType value{};
        insert(key, value);
        it = lookup(key);
        if (it == nullptr) {
            it = lookup(key);
        }
    }
    return (*it)->second;
}

template<class KeyType, class ValueType, class Hash>
const ValueType &HashMap<KeyType, ValueType, Hash>::at(const KeyType &key) const {
    auto it = lookup_const(key);
    if (it == nullptr) {
        throw std::out_of_range("Key is not in HashMap");
    }
    return (*it)->second;
}


/**********************************************************************
 *                         begin(), end()                             *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::const_iterator HashMap<KeyType, ValueType, Hash>::end() const {
    return const_iterator{positions_.end()};
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::const_iterator HashMap<KeyType, ValueType, Hash>::begin() const {
    return const_iterator{positions_.begin()};
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::iterator HashMap<KeyType, ValueType, Hash>::begin() {
    return iterator{positions_.begin()};
}

template<class KeyType, class ValueType, class Hash>
typename HashMap<KeyType, ValueType, Hash>::iterator HashMap<KeyType, ValueType, Hash>::end() {
    return iterator{positions_.end()};
}

/**********************************************************************
 *                          HashFunction                              *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::HashFunction::randomize() {
    y = rng();
    k = rng();
    b = rng();
}

template<class KeyType, class ValueType, class Hash>
size_t HashMap<KeyType, ValueType, Hash>::HashFunction::operator()(size_t x, size_t p) const {
    x ^= y;
    x %= p;
    x *= k % p;
    x %= p;
    x += b % p;
    if (x >= p) {
        x -= p;
    }
    return x;
}

/**********************************************************************
 *                       Helper functions                             *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
size_t HashMap<KeyType, ValueType, Hash>::max_sum_s() const {
//    return 32 * m_ * m_ / S(m_) + 4 * m_;
    return 25 * m_;
}

template<class KeyType, class ValueType, class Hash>
size_t HashMap<KeyType, ValueType, Hash>::S(size_t m) {
    return static_cast<size_t>(0.7 * m);
}

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::fix_j(
        std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>> &l_j,
        std::vector<size_t> &hashes, size_t j) {
    std::vector<char> has(t_[j].size());
    while (true) {
        hash_local_[j].randomize();
        bool bad = false;
        for (size_t hash: hashes) {
            if (has[hash_local_[j](hash, t_[j].size())]) {
                bad = true;
                break;
            }
            has[hash_local_[j](hash, t_[j].size())] = true;
        }
        if (!bad) {
            for (size_t i = 0; i < l_j.size(); ++i) {
                auto h = hash_local_[j](hashes[i], t_[j].size());
                t_[j][h].values.swap(l_j[i]);
                auto it = t_[j][h].values.begin();
                while (it != t_[j][h].values.end()) {
                    it->second->iter_ = it;
                    ++it;
                }
            }
            break;
        }
        has.assign(has.size(), 0);
    }
}

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::rehash_all(const std::vector<std::pair<const KeyType, ValueType>> &inserted) {
    std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>> l;
    std::vector<size_t> hashes;
    for (size_t i = 0; i < t_.size(); ++i) {
        for (size_t j = 0; j < t_[i].size(); ++j) {
            if (!t_[i][j].values.empty()) {
                l.push_back({});
                t_[i][j].values.swap(l.back());
                hashes.push_back(hasher_(l.back().begin()->first.first));
            }
        }
    }
    for (auto &elem: inserted) {
        bool good = false;
        auto cur_hash = hasher_(elem.first);
        for (size_t i = 0; i < l.size(); ++i) {
            if (cur_hash == hashes[i]) {
                positions_.push_back({{}});
                l[i].push_back({elem, --positions_.end()});
                good = true;
            }
        }
        if (!good) {
            positions_.push_back({{}});
            l.push_back({{elem, --positions_.end()}});
            hashes.push_back(hasher_(l.back().begin()->first.first));
        }
    }
    count_ = size_;
    m_ = (1 + HASH_MAP_C) * std::max(count_, 4UL);
    while (true) {
        hash_global_.randomize();
        t_.clear();
        t_.resize(get_nearest_prime(S(m_)));
        hash_local_.resize(t_.size());
        m_local_.assign(t_.size(), 0);
        s_local_.assign(t_.size(), 0);
        b_.assign(t_.size(), 0);
        for (size_t hash: hashes) {
            ++b_[hash_global_(hash, t_.size())];
        }
        sum_s_ = 0;
        for (size_t i = 0; i < t_.size(); ++i) {
            if (b_[i] == 0) {
                m_local_[i] = 1;
                s_local_[i] = 1;
            } else {
                m_local_[i] = std::max(2UL, static_cast<size_t>(2 * b_[i]));
                s_local_[i] = get_nearest_prime(3 * m_local_[i] * (m_local_[i] - 1));
            }
            sum_s_ += s_local_[i];

        }
        if (sum_s_ <= max_sum_s()) {
            std::vector<std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>>> l_j(
                    t_.size());
            std::vector<std::vector<size_t>> hashes_j(t_.size());
            for (size_t i = 0; i < l.size(); ++i) {
                size_t j = hash_global_(hashes[i], t_.size());
                l_j[j].push_back(l[i]);
                hashes_j[j].push_back(hashes[i]);
            }
            for (size_t i = 0; i < t_.size(); ++i) {
                t_[i].resize(s_local_[i]);
                fix_j(l_j[i], hashes_j[i], i);
            }
            break;
        }
    }
}

/**********************************************************************
 *                        lookup helpers                              *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
std::unique_ptr<typename HashMap<KeyType, ValueType, Hash>::const_iterator>
HashMap<KeyType, ValueType, Hash>::lookup_const(const KeyType &key) const {
    if (t_.empty()) {
        return nullptr;
    }
    auto h = hasher_(key);
    size_t global_index = hash_global_(h, t_.size());
    size_t local_index = hash_local_[global_index](h, t_[global_index].size());
    bool good = false;
    for (size_t i = 0; i < MAX_OFFSET; ++i) {
        size_t index = (local_index + i) % t_[global_index].size();
        if (!t_[global_index][index].values.empty() &&
            hasher_(t_[global_index][index].values.begin()->first.first) == h) {
            local_index = index;
            good = true;
            break;
        }
    }
    if (!good) {
        return nullptr;
    }
    typename std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>::const_iterator iter = t_[global_index][local_index].values.begin();
    while (iter != t_[global_index][local_index].values.end() && !(iter->first.first == key)) { // NOLINT
        ++iter;
    }
    if (iter == t_[global_index][local_index].values.end()) {
        return nullptr;
    }
    return std::make_unique<const_iterator>(iter->second);
}

template<class KeyType, class ValueType, class Hash>
std::unique_ptr<typename HashMap<KeyType, ValueType, Hash>::iterator>
HashMap<KeyType, ValueType, Hash>::lookup(const KeyType &key) {
    if (t_.empty()) {
        return nullptr;
    }
    auto h = hasher_(key);
    size_t global_index = hash_global_(h, t_.size());
    size_t local_index = hash_local_[global_index](h, t_[global_index].size());
    bool good = false;
    for (size_t i = 0; i < MAX_OFFSET; ++i) {
        size_t index = (local_index + i) % t_[global_index].size();
        if (!t_[global_index][index].values.empty() &&
            hasher_(t_[global_index][index].values.begin()->first.first) == h) {
            local_index = index;
            good = true;
            break;
        }
    }
    if (!good) {
        return nullptr;
    }
    typename std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>::const_iterator iter = t_[global_index][local_index].values.begin();
    while (iter != t_[global_index][local_index].values.end() && !(iter->first.first == key)) { // NOLINT
        ++iter;
    }
    if (iter == t_[global_index][local_index].values.end()) {
        return nullptr;
    }
    return std::make_unique<iterator>(iter->second);
}

/**********************************************************************
 *                  Some other required functions                     *
 **********************************************************************/


template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::contains(const KeyType &key) {
    return static_cast<bool>(lookup_const(key));
}

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::insert(const std::pair<KeyType, ValueType> &key_value) {
    insert(key_value.first, key_value.second);
}

template<class KeyType, class ValueType, class Hash>
const Hash &HashMap<KeyType, ValueType, Hash>::hash_function() const {
    return hasher_;
}

template<class KeyType, class ValueType, class Hash>
size_t HashMap<KeyType, ValueType, Hash>::size() const {
    return size_;
}


template<class KeyType, class ValueType, class Hash>
bool HashMap<KeyType, ValueType, Hash>::empty() const {
    return size_ == 0;
}

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::clear() {
    HashMap<KeyType, ValueType, Hash> tmp(hasher_);
    swap(tmp);
}

/**********************************************************************
 *                             insert                                 *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::insert(const KeyType &key, const ValueType &value) {
    size_t key_hash = hasher_(key);
    if (contains(key)) {
        return;
    }
    ++size_;
    if (t_.empty()) {
        rehash_all({{key, value}});
        return;
    }
    size_t j = hash_global_(hasher_(key), t_.size());
    auto local = hash_local_[j](key_hash, t_[j].size());
    if (!t_.empty()) {
        if (!t_[j][local].values.empty() && hasher_(t_[j][local].values.begin()->first.first) == key_hash) {
            t_[j][local].values.push_back({{key, value},
                                           {}});
            positions_.push_back({--t_[j][local].values.end()});
            t_[j][local].values.back().second = --positions_.end();
            return;
        }
    }
    {
        // kinda like open addressing
        for (size_t i = 0; i < MAX_OFFSET; ++i) {
            size_t pos = (local + i) % t_[j].size();
            if (!t_[j][pos].values.empty() && hasher_(t_[j][pos].values.begin()->first.first) == key_hash) {
                t_[j][pos].values.push_back({{key, value},
                                             {}});
                positions_.push_back({--t_[j][pos].values.end()});
                t_[j][pos].values.back().second = --positions_.end();
                return;
            }
        }
        for (size_t i = 0; i < MAX_OFFSET; ++i) {
            size_t pos = (local + i) % t_[j].size();
            if (t_[j][pos].values.empty()) {
                t_[j][pos].values.push_back({{key, value},
                                             {}});
                positions_.push_back({--t_[j][pos].values.end()});
                t_[j][pos].values.back().second = --positions_.end();
                b_[j]++;
                return;
            }
        }
    }
    if (size_ > m_ * 2.6) {
        rehash_all({{key, value}});
        return;
    }
    ++count_;
    if (count_ > m_) {
        rehash_all({{key, value}});
    } else {
        ++b_[j];
        if (b_[j] <= m_local_[j]) {
            if (t_[j][local].values.empty()) {
                t_[j][local].values.push_back({{key, value},
                                               {}});
                positions_.push_back({--t_[j][local].values.end()});
                t_[j][local].values.back().second = --positions_.end();
            } else {
                std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>> l_j;
                std::vector<size_t> hashes;
                l_j.reserve(t_[j].size());
                hashes.reserve(t_[j].size());
                for (auto &el: t_[j]) {
                    if (!el.values.empty()) {
                        l_j.push_back({});
                        el.values.swap(l_j.back());
                        hashes.push_back(hasher_(l_j.back().begin()->first.first));
                    }
                }
                positions_.push_back({{}});
                l_j.push_back({{{key, value}, --positions_.end()}});
                hashes.push_back(key_hash);
                b_[j] = hashes.size();
                fix_j(l_j, hashes, j);
            }
        } else {
            m_local_[j] = 2 * std::max(1UL, m_local_[j]);
            sum_s_ -= s_local_[j];
            s_local_[j] = get_nearest_prime(2 * m_local_[j] * (m_local_[j] - 1));
            sum_s_ += s_local_[j];
            if (sum_s_ <= max_sum_s()) {
                std::vector<std::list<std::pair<std::pair<const KeyType, ValueType>, typename std::list<Position>::iterator>>> l_j;
                std::vector<size_t> hashes;
                l_j.reserve(t_[j].size());
                hashes.reserve(t_[j].size());
                for (auto &el: t_[j]) {
                    if (!el.values.empty()) {
                        l_j.push_back({});
                        el.values.swap(l_j.back());
                        hashes.push_back(hasher_(l_j.back().begin()->first.first));
                    }
                }
                positions_.push_back({{}});
                l_j.push_back({{{key, value}, --positions_.end()}});
                hashes.push_back(key_hash);
                t_[j].resize(s_local_[j]);
                b_[j] = hashes.size();
                fix_j(l_j, hashes, j);
            } else {
                rehash_all({{key, value}});
            }
        }
    }
}

/**********************************************************************
 *                              erase                                 *
 **********************************************************************/

template<class KeyType, class ValueType, class Hash>
void HashMap<KeyType, ValueType, Hash>::erase(const KeyType &key) {
    if (t_.empty()) {
        return;
    }
//    ++count_;
    auto h = hasher_(key);
    size_t global_index = hash_global_(h, t_.size());
    size_t local_index = hash_local_[global_index](h, t_[global_index].size());
    for (size_t i = 0; i < MAX_OFFSET; ++i) {
        size_t index = (local_index + i) % t_[global_index].size();
        if (!t_[global_index][index].values.empty() &&
            hasher_(t_[global_index][index].values.begin()->first.first) == h) {
            local_index = index;
            break;
        }
    }
    auto it = t_[global_index][local_index].values.begin();
    while (it != t_[global_index][local_index].values.end() && !(it->first.first == key)) { // NOLINT
        ++it;
    }
    if (it == t_[global_index][local_index].values.end()) {
        return;
    }
    --size_;
    positions_.erase(it->second);
    t_[global_index][local_index].values.erase(it);
    if (count_ > m_) {
//        rehash_all({});
    }
}
