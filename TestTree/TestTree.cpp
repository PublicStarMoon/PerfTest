// Statistics.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <chrono>
#include <map>
#include <set>
#include <vector>
#include <array>
#include <tuple>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <atomic>

#include <cmath>

#include "avlmini.h"
#include "test_avl.h"
#include "test/test_linux_rb.h"


//---------------------------------------------------------------------
// random 
//---------------------------------------------------------------------

// (Mean, Variance)
using Statistics = std::tuple<double, double>;
constexpr size_t NODE_COUNT = 1000000;
constexpr size_t ITERATIONS = 1000;

std::mutex io_mutex;

template<typename RandomAccessIterator>
Statistics
CalculateStatistic(RandomAccessIterator begin, RandomAccessIterator end)
{
	auto size = std::distance(begin, end);

	double sum = 0;
	for (auto i = begin; i != end; ++i)
	{
		sum += *i;
	}

	double mean{ 0 }, variance{ 0 };

	mean = sum / size;

	sum = 0;
	for (auto i = begin; i != end; ++i)
	{
		sum += (*i - mean) * (*i - mean);
	}
	variance = std::sqrt(sum / size);

	return std::make_tuple(mean, variance);
}

template<typename RandomAccessContainer>
void SampleRbTree(RandomAccessContainer& insert, RandomAccessContainer& find, RandomAccessContainer& erase)
{
	{
		std::lock_guard<std::mutex> lock(io_mutex);
		std::cout << "SampleRbTree" << std::endl;
	}

	struct rb_root rb_root;
	rb_root.rb_node = nullptr;

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();

	std::vector<int> data(NODE_COUNT, 0);
	std::iota(data.begin(), data.end(), 0);

	std::vector<struct rb_node*> rb_pool(NODE_COUNT, nullptr);
	for (size_t i = 0; i != NODE_COUNT; ++i) {
		rb_pool[i] = reinterpret_cast<struct rb_node*>(avl_node_new(static_cast<int>(i)));
	}

	for (size_t it = 0; it != ITERATIONS; ++it)
	{
		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample insert
		for (auto i = data.begin(); i != data.end(); ++i) {
			struct rb_node* dup = nullptr;
			rb_node_add(&rb_root, rb_pool[*i], rb_node_compare, dup);
		}
		end = std::chrono::system_clock::now();
		insert[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();

		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample find
		for (auto i = data.begin(); i != data.end(); ++i) {
			struct rb_node* res = nullptr;
			struct RbNode dummy;
			dummy.key = *i;
			rb_node_find(&rb_root, &dummy.node, rb_node_compare, res);

			struct MyNode* result = AVL_ENTRY(res, struct MyNode, node);
			//assert(result);
			//assert(result->key == *i);
			if (!result || result->key != *i) {
				std::cout << (res == rb_pool[*i]) << std::endl;
			}
		}
		end = std::chrono::system_clock::now();
		find[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();

		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample erase
		for (auto i = data.begin(); i != data.end(); ++i) {
			//struct rb_node* res = nullptr;
			//struct RbNode dummy;
			//dummy.key = *i;
			//rb_node_find(&rb_root, &dummy.node, rb_node_compare, res);
			rb_erase(rb_pool[*i], &rb_root);
		}
		end = std::chrono::system_clock::now();
		erase[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();

		std::cout << it << std::endl;
	}

	auto insert_statistics = CalculateStatistic(insert.begin(), insert.end());
	auto find_statistics = CalculateStatistic(find.begin(), find.end());
	auto erase_statistics = CalculateStatistic(erase.begin(), erase.end());

	{
		std::lock_guard<std::mutex> lock(io_mutex);
		std::cout << "Insert (mean, variance) ns: "
			<< std::get<0>(insert_statistics) << " , " << std::get<1>(insert_statistics) << std::endl;
		std::cout << "Find (mean, variance) ns: "
			<< std::get<0>(find_statistics) << " , " << std::get<1>(find_statistics) << std::endl;
		std::cout << "Erase (mean, variance) ns: "
			<< std::get<0>(erase_statistics) << " , " << std::get<1>(erase_statistics) << std::endl;
	}

	for (size_t i = 0; i != NODE_COUNT; ++i) {
		free(rb_pool[i]);
	}
}

template<typename RandomAccessContainer>
void SampleAvlTree(RandomAccessContainer& insert, RandomAccessContainer& find, RandomAccessContainer& erase)
{
	{
		std::lock_guard<std::mutex> lock(io_mutex);
		std::cout << "SampleAvlTree" << std::endl;
	}

	struct avl_root avl_root;
	avl_root.node = nullptr;

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();

	std::vector<int> data(NODE_COUNT, 0);
	std::iota(data.begin(), data.end(), 0);

	std::vector<struct avl_node*> avl_pool(NODE_COUNT, nullptr);
	for (size_t i = 0; i != NODE_COUNT; ++i) {
		avl_pool[i] = reinterpret_cast<struct avl_node*>(avl_node_new(static_cast<int>(i)));
	}

	for (size_t it = 0; it != ITERATIONS; ++it)
	{
		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample insert
		for (auto i = data.begin(); i != data.end(); ++i) {
			struct avl_node* dup = nullptr;
			avl_node_add(&avl_root, avl_pool[*i], avl_node_compare, dup);
		}
		end = std::chrono::system_clock::now();
		insert[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();

		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample find
		for (auto i = data.begin(); i != data.end(); ++i) {
			struct avl_node* res = nullptr;
			struct MyNode dummy;
			dummy.key = *i;
			avl_node_find(&avl_root, &dummy.node, avl_node_compare, res);

			struct MyNode* result = AVL_ENTRY(res, struct MyNode, node);
			if (!result || result->key != *i) {
				std::cout << (res == avl_pool[*i]) << std::endl;
			}
		}
		end = std::chrono::system_clock::now();
		find[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();

		std::random_shuffle(data.begin(), data.end());

		start = std::chrono::system_clock::now();
		// sample erase
		for (auto i = data.begin(); i != data.end(); ++i) {
			avl_node_erase(avl_pool[*i], &avl_root);
		}
		end = std::chrono::system_clock::now();
		erase[it] = std::chrono::duration_cast<std::chrono::nanoseconds>((end - start)).count();
	}

	auto insert_statistics = CalculateStatistic(insert.begin(), insert.end());
	auto find_statistics = CalculateStatistic(find.begin(), find.end());
	auto erase_statistics = CalculateStatistic(erase.begin(), erase.end());

	{
		std::lock_guard<std::mutex> lock(io_mutex);
		std::cout << "Insert (mean, variance) ns: "
			<< std::get<0>(insert_statistics) << " , " << std::get<1>(insert_statistics) << std::endl;
		std::cout << "Find (mean, variance) ns: "
			<< std::get<0>(find_statistics) << " , " << std::get<1>(find_statistics) << std::endl;
		std::cout << "Erase (mean, variance) ns: "
			<< std::get<0>(erase_statistics) << " , " << std::get<1>(erase_statistics) << std::endl;
	}

	for (size_t i = 0; i != NODE_COUNT; ++i) {
		free(avl_pool[i]);
	}
}

template<typename RandomAccessContainer>
std::vector<std::thread>
createAvlThread(
	RandomAccessContainer& insert,
	RandomAccessContainer& find,
	RandomAccessContainer& erase,
	size_t thread_num)
{
	for (size_t i = 0; i < thread_num; ++i)
	{
		auto workerInsertPtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);
		auto workerFindPtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);
		auto workerErasePtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);

		insert.emplace_back(workerInsertPtr);
		find.emplace_back(workerFindPtr);
		erase.emplace_back(workerErasePtr);
	}

	std::vector<std::thread> pool;
	for (size_t i = 0; i < thread_num; ++i)
	{
		pool.emplace_back(std::thread(
			SampleAvlTree<std::vector<long long>>,
			std::reference_wrapper<std::vector<long long>>(*(insert[i])),
			std::reference_wrapper<std::vector<long long>>(*(find[i])),
			std::reference_wrapper<std::vector<long long>>(*(erase[i]))
		));
	}
	return pool;
}

template<typename RandomAccessContainer>
std::vector<std::thread>
createRbThread(
	RandomAccessContainer& insert,
	RandomAccessContainer& find,
	RandomAccessContainer& erase,
	size_t thread_num)
{
	for (size_t i = 0; i < thread_num; ++i)
	{
		auto workerInsertPtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);
		auto workerFindPtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);
		auto workerErasePtr = std::make_shared<std::vector<long long>>(ITERATIONS, 0);

		insert.emplace_back(workerInsertPtr);
		find.emplace_back(workerFindPtr);
		erase.emplace_back(workerErasePtr);
	}

	std::vector<std::thread> pool;
	for (size_t i = 0; i < thread_num; ++i)
	{
		pool.emplace_back(std::thread(
			SampleRbTree<std::vector<long long>>,
			std::reference_wrapper<std::vector<long long>>(*(insert[i])),
			std::reference_wrapper<std::vector<long long>>(*(find[i])),
			std::reference_wrapper<std::vector<long long>>(*(erase[i]))
		));
	}
	return pool;
}

template<typename ForwardIterator>
void SaveToCsvMultiThread(
	std::string csv_path,
	ForwardIterator insertBegin,
	ForwardIterator insertEnd,
	ForwardIterator findBegin,
	ForwardIterator findEnd,
	ForwardIterator eraseBegin,
	ForwardIterator eraseEnd)
{
	std::ofstream fout(csv_path, std::ios_base::out | std::ios_base::trunc);

	fout << "insert,find,erase" << std::endl;

	using FwdIterator = typename ForwardIterator::value_type::element_type::iterator;

	ForwardIterator i, j, k;
	for (i = insertBegin, j = findBegin, k = eraseBegin;
		i != insertEnd && j != findEnd && k != eraseEnd;
		++i, ++j, ++k)
	{
		FwdIterator ii, jj, kk;
		for (ii = (*i)->begin(), jj = (*j)->begin(), kk = (*k)->begin();
			ii != (*i)->end() && jj != (*j)->end() && kk != (*k)->end();
			++ii, ++jj, ++kk)
		{
			fout << *ii << "," << *jj << "," << *kk << std::endl;
		}
	}

	fout.close();
}

int main()
{
	std::vector<std::shared_ptr<std::vector<long long>>> avl_insert;
	std::vector<std::shared_ptr<std::vector<long long>>> avl_find;
	std::vector<std::shared_ptr<std::vector<long long>>> avl_erase;

	std::vector<std::shared_ptr<std::vector<long long>>> rb_insert;
	std::vector<std::shared_ptr<std::vector<long long>>> rb_find;
	std::vector<std::shared_ptr<std::vector<long long>>> rb_erase;

	std::vector<std::thread> avl_pool = std::move(createAvlThread(avl_insert, avl_find, avl_erase, 10));
	std::vector<std::thread> rb_pool = std::move(createRbThread(rb_insert, rb_find, rb_erase, 10));

	for (auto&& t : rb_pool) t.join();
	for (auto&& t : avl_pool) t.join();

	SaveToCsvMultiThread("avltree.csv", 
		avl_insert.begin(), avl_insert.end(), avl_find.begin(), avl_find.end(), avl_erase.begin(), avl_erase.end());
	SaveToCsvMultiThread("rbtree.csv",
		rb_insert.begin(), rb_insert.end(), rb_find.begin(), rb_find.end(), rb_erase.begin(), rb_erase.end());
}