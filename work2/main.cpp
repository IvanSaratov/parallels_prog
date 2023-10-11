#include <iostream>
#include <time.h>
#include <omp.h>
#include <Windows.h>
#include <cstdio>

using namespace std;

constexpr auto PI = 3.1415926535897932384626433832795;

long double func(const double x, const double y, const double a1, const double b1, const double a2, const double b2);

void integralCells(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res);

long double experiment(double* result);

int main() {
	SetConsoleOutputCP(1251);   // Устанавливаем кодировку для windows


	long long i; // переменная цикла
	long double time; // время проведенного эксперимента
	double result; // значение вычисленного интеграла
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

long double func(const double x, const double y, const double a1, const double b1, const double a2, const double b2) {
	return (exp(sin(PI * x) * cos(PI * y)) + 1.0) / ((b1 - a1) * (b2 - a2));
}

void integralCells(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* result) {
	long long i, j;
	double sum = 0.0; //Начальное значение интеграла
	double x0 = a1, y0 = a2, x1, y1;
	long n = (long)((b1 - a1) / h1); // dx
	long m = (long)((b2 - a2) / h2); // dy

	#pragma omp parallel for private(x1, y1, x0, y0) reduction(+: sum)
	for (i = 1; i < n + 1; i++) {
		x1 = x0 + h1;
		for (j = 1; j < m + 1; j++) {
			y1 = y0 + h2;
			sum += func((x0 + x1) / 2.0, (y0 + y1) / 2.0, a1, b1, a2, b2);
			y0 = y1;
		}
		x0 = x1;
	}
	sum *= h1 * h2;
	*result = sum;
}

long double experiment(double* result) {
	double stime, ftime; // Время начала и конца расчета
	double a1 = 0.0;
	double a2 = 0.0;
	double b1 = 16.0;
	double b2 = 16.0;
	double h1 = 0.01; // Шаг интегрирования
	double h2 = 0.01; // Шаг интегрирования

	stime = clock();
	integralCells(a1, b1, a2, b2, h1, h2, result); // Вычисляем
	ftime = clock();
	return (ftime - stime) / CLOCKS_PER_SEC;
}