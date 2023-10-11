#include <iostream>
#include <time.h>
#include <omp.h>
#include <Windows.h>
#include <cstdio>

using namespace std;

// Определяем PI
constexpr auto PI = 3.1415926535897932384626433832795;

// Заменяем бесконечность
const long INF = 1000000;

double func(double x);

void integralSimple(const long double a, const long double b, const long double h, long double* result);

void integralSimpsona(const long double a, const long double b, const long long n, long double* result);

double experiment(long double* res);

int main()
{
	SetConsoleOutputCP(1251);   // Устанавливаем кодировку для windows


	long long i; // переменная цикла
	long double time; // время проведенного эксперимента
	long double result; // значение вычисленного интеграла
	long double min_time; // минимальное время работы реализации алгоритма
	long double max_time; // максимальное время работы реализации алгоритма
	long double avg_time; // среднее время работы реализации алгоритма
	long long numbExp = 5; // количество запусков программы

	// Первый запуск
	min_time = max_time = avg_time = experiment(&result);
	// Остальные запуски
	for (i = 0; i < numbExp - 1; i++)
	{
		time = experiment(&result);
		avg_time += time;
		if (max_time < time) max_time = time;
		if (min_time > time) min_time = time;
	}

	// Вывод работы эксперемента
	// Минимальное - среднее - максимальное
	cout.precision(8);
	cout << "Время выполнения: " << min_time << "; " << avg_time / numbExp << "; " << max_time << endl;
	cout << "Интегральное значение : " << result << endl;
	cout << "Ответ который должен быть: " << PI / 2.0 << endl;

	return 0;
}

// Наша заданная функция по варианту
double func(double x)
{
	return x / (1.0 + x * x);
}

void integralSimple(const long double a, const long double b, const long double h, long double* result)
{
	long long i, j = 0;
	double sum = 0.0; // Локальная переменная для подсчета интеграла
	double x = 0.0;  // Кордината точки сетки
	// Определяем количество точек сетки интегрирования

	j = (long long)((b - a) / h);

	#pragma omp parallel for private(x) reduction(+: sum) 
	for (i = 1; i <= j; i++)
	{
		x = a + i * h + h / 2.0;
		sum += func(x) * h;
	}
	*result = sum;
}

void integralSimpsona(const long double a, const long double b, const long long n, long double* result) {
	long double h = 0.0;
	double k = 0.0, s = 0.0;
	double x = 0.0;

	h = (b - a) / (2 * n);
	#pragma omp parallel for private(x) reduction(+: k)
	for (long long i = 1; i <= 2 * n - 1; i++)
	{
		if (i % 2 == 0)
			k += 2 * func(a + i * h);
		else 
			k += 4 * func(a + i * h);
	}
	s = ((func(a) + func(b) + k)) * h / 3;

	*result = s;
}


double experiment(long double* result)
{
	long double stime, ftime; // Время начала и конца расчета
	long double a = 0.0; // Левая граница интегрирования
	long double b = INF; // Правая граница интегрирования
	long double h = 1e9; // Шаг интегрирования

	stime = clock();
	//integralSimple(a, b, h, result); // Вызов функции интегрирования
	integralSimpsona(a, b, h, result); // Вызов метода Симпсона
	ftime = clock();
	return (ftime - stime) / CLOCKS_PER_SEC;
}