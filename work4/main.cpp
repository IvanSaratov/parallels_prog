#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <time.h>
#include <omp.h>
#include <Windows.h>
#include <cstdio>

using namespace std;

constexpr auto PI = (3.14159265358979323846);

// Выделяем память и инициализируем значения
void ProcessInitialization(complex<double>*& inputSignal, complex<double>*& outputSignal, int& size);
// Функция заполнения случайными значениями
void RandomDataInitialization(complex<double>* mas, int size);

void BitReversing(complex<double>* inputSignal, complex<double>* outputSignal, int size);


void BitReversingParallel(complex<double>* inputSignal, complex<double>* outputSignal, int size);
__inline void Butterfly(complex<double>* signal, complex<double> u, int offset, int butterflySize);

// Запус алгоритма БФТ последовательно
void SingleFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size);
// Логика алгоритма
void SingleFFTCalculation(complex<double>* signal, int size);

// Запуск параллельного алгоритма
void ParallelFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size);
// Логика работы параллельного алгоритма
void ParallelFFTCalculation(complex<double>* signal, int size);


// Функция выода нашего сигнала на экран
void PrintSignal(complex<double>* signal, int size);
// Функция запуска сравнения правильности вводных данных
void TestResult(complex<double>* inputSignal, complex<double>* outputSignal, int size);

// Очистка памяти
void ProcessTermination(complex<double>*& inputSignal, complex<double>*& outputSignal);

int main()
{
	SetConsoleOutputCP(1251);   // Устанавливаем кодировку для windows

	complex<double>* inputSignal = NULL;
	complex<double>* outputSignal = NULL;

	int size = 0;
	const int repeatCount = 16;
	double startTime;
	double duration;
	double minDuration = DBL_MAX;

	cout << "Быстрое преобразование Фурье" << endl;
	// Инициализируем наши сигналы и заполняем их
	ProcessInitialization(inputSignal, outputSignal, size);

	for (int i = 0; i < repeatCount; i++)
	{
		startTime = clock();

		// Запускаем вычисление в одном потоке
		//SingleFFT(inputSignal, outputSignal, size);
		// Запускаем вычисление в несколько потоков
		ParallelFFT(inputSignal, outputSignal, size);

		// Замеряем время
		duration = (clock() - startTime) / CLOCKS_PER_SEC;
		if (duration < minDuration)
			minDuration = duration;
	}

	cout << setprecision(6);
	cout << endl << "Время выполнения: " << minDuration << " секунд" << endl;

	// Запускам наши тесты
	TestResult(inputSignal, outputSignal, size);

	// Очищаем память и завершаем программу
	ProcessTermination(inputSignal, outputSignal);
	return 0;
}


void ProcessInitialization(complex<double>*& inputSignal, complex<double>*& outputSignal, int& size) {
	// Устанавливаем длину сигнала (Например 8)
	do
	{
		cout << "Введите длину сигнала: ";
		cin >> size;
		if (size < 4)
			cout << "Длина сигнала должны быть больше 4" << endl;
		else
		{
			int tmpSize = size;
			while (tmpSize != 1)
			{
				if (tmpSize % 2 != 0)
				{
					cout << "Длина сигнала должна быть кратна 2" << endl;
					size = -1;
					break;
				}
				tmpSize /= 2;
			}
		}
	} while (size < 4);

	cout << "Введенная длина сигнала = " << size << endl;
	inputSignal = new complex<double>[size];
	outputSignal = new complex<double>[size];

	// Заполняем случайными значениями
	RandomDataInitialization(inputSignal, size);
}

void RandomDataInitialization(complex<double>* mas, int size)
{
	srand(unsigned(clock()));
	for (int i = 0; i < size; i++)
		mas[i] = complex<double>(rand() / 1000.0, rand() / 1000.0);
}

void PrintSignal(complex<double>* signal, int size) {
	cout << "Сигнал: " << endl;
	for (int i = 0; i < size; i++)
		cout << signal[i] << endl;
}

void SingleFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	BitReversing(inputSignal, outputSignal, size);
	SingleFFTCalculation(outputSignal, size);
}

void SingleFFTCalculation(complex<double>* signal, int size) {
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);
	for (int p = 0; p < m; p++) {
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;
		double coeff = PI / butterflySize;
		for (int i = 0; i < size / butterflyOffset; i++)
			for (int j = 0; j < butterflySize; j++)
				Butterfly(signal, complex<double>(cos(-j * coeff), sin(-j * coeff)), j + i * butterflyOffset, butterflySize);
	}
}

void ParallelFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	BitReversing(inputSignal, outputSignal, size);
	ParallelFFTCalculation(outputSignal, size);
}

void ParallelFFTCalculation(complex<double>* signal, int size) {
	int m = 0, i = 0, j = 0;
	int butterflyOffset = 0;
	int butterflySize = 0;

	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);

	for (int p = 0; p < m; p++)
	{
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;
		double coeff = PI / butterflySize;

#pragma omp parallel for private(i, j)
		for (i = 0; i < size / butterflyOffset; i++)
			for (j = 0; j < butterflySize; j++)
				Butterfly(signal, complex<double>(cos(-j * coeff), sin(-j * coeff)), j + i * butterflyOffset, butterflySize);
	}
}

void BitReversing(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	int j = 0, i = 0;

	while (i < size)
	{
		if (j > i)
		{
			outputSignal[i] = inputSignal[j];
			outputSignal[j] = inputSignal[i];
		}
		else if (j == i)
			outputSignal[i] = inputSignal[i];
		int m = size >> 1;
		while ((m >= 1) && (j >= m)) {
			j -= m;
			m = m >> 1;
		}
		j += m;
		i++;
	}
}

void BitReversingParallel(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	int j = 0, i = 0;

#pragma omp parallel for ordered
	for (i = 0; i < size; i++)
	{
		if (j > i)
		{
			outputSignal[i] = inputSignal[j];
			outputSignal[j] = inputSignal[i];
		}
		else
			if (j == i)
				outputSignal[i] = inputSignal[i];
		int m = size >> 1;
		while ((m >= 1) && (j >= m))
		{
			j -= m;
			m = m >> 1;
		}
#pragma omp ordered
		j += m;
	}
}

__inline void Butterfly(complex<double>* signal, complex<double> u, int offset, int butterflySize) {
	complex<double> tem = signal[offset + butterflySize] * u;
	signal[offset + butterflySize] = signal[offset] - tem;
	signal[offset] += tem;
}

void TestResult(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	// Буффер для временного хранения сигналов
	complex<double>* testSerialSignal;
	double Accuracy = 1.e-6;
	int equal = 0;
	int i;
	testSerialSignal = new complex<double>[size];
	// Запускам проверку в одном потоке
	SingleFFT(inputSignal, testSerialSignal, size);

	for (i = 0; i < size; i++)
	{
		// Сравниваем параллельное выполнение и последовательное
		if (abs(outputSignal[i] - testSerialSignal[i]) >= Accuracy)
			equal = 1;
	}
	if (equal == 1)
		cout << "Результаты НЕ совпали" << endl;
	else
		cout << "Результаты совпали" << endl;

	// Освобождаем место из под буффера
	delete[] testSerialSignal;
}

void ProcessTermination(complex<double>*& inputSignal, complex<double>*& outputSignal) {
	delete[] inputSignal;
	inputSignal = NULL;
	delete[] outputSignal;
	outputSignal = NULL;
}