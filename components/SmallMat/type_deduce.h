#pragma once

namespace matrix{

	// type guess
	template<typename t1, typename t2>class TypeDeduce{};

	template<> class TypeDeduce<int, int>{public:typedef int type;};
	template<> class TypeDeduce<int, double>{public:typedef double type;};
	template<> class TypeDeduce<int, std::complex<double>>{
		public:typedef std::complex<double> type;};

	template<> class TypeDeduce<double, int>{public:typedef double type;};
	template<> class TypeDeduce<double, double>{public:typedef double type;};
	template<> class TypeDeduce<double, std::complex<double>>{
		public:typedef std::complex<double> type;};

	template<> class TypeDeduce<std::complex<double>, int>{
		public:typedef std::complex<double> type;};
	template<> class TypeDeduce<std::complex<double>, double>{
		public:typedef std::complex<double> type;};
	template<> class TypeDeduce<std::complex<double>, std::complex<double>>{
		public:typedef std::complex<double> type;};

};