#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <Windows.h>
#include <cstdio>

using namespace std;

typedef struct {
	int PivotRow;
	double MaxValue;
} TThreadPivotRow;

int* pPivotPos;  // The number of pivot rows selected at the iterations
int* pPivotIter; // The iterations, at which the rows were pivots


// Выделяем память для нашей матрицы и заполняем её значениями
void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int& Size);
// Функция заполнения предзаготовленными параметрами
void DataInitializationPreset(double* pMatrix, double* pVector, int Size);
// Функция для вывода матрицы на экран
void PrintMatrix(double* pMatrix, int RowCount, int ColCount);
// Функция для вывода вектора на экран
void PrintVector(double* pVector, int Size);


// Функции запуска алгоритма Гаусса
void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size);
int ParallelFindPivotRow(double* pMatrix, int Size, int Iter);
void ParallelColumnElimination(double* pMatrix, double* pVector, int Pivot, int Iter, int Size);
void ParallelGaussianElimination(double* pMatrix, double* pVector, int Size);
void ParallelBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size);

// Функция для проверки парвильность вычисления методом Гаусса
void TestResult(double* pMatrix, double* pVector, double* pResult, int Size);
// Функция очистки памяти
void ProcessTermination(double* pMatrix, double* pVector, double* pResult);

int main() {
	SetConsoleOutputCP(1251);   // Устанавливаем кодировку для windows

	double* pMatrix; // Наша матрица
	double* pVector; // Правая сторона матрицы
	double* pResult; // Результат
	int Size = 5; // Размер матрицы. Заранее указываем
	double start, finish, duration;

	cout << "Параллельный алгоритм Гаусса для линейных систем" << endl;
	cout << "Загружаем пресет значений" << endl;

	ProcessInitialization(pMatrix, pVector, pResult, Size);

	// The matrix and the vector output
	cout << endl << "Наша матрица: " << endl;
	PrintMatrix(pMatrix, Size, Size);
	cout << endl << "Права сторона нашей матрицы: " << endl;
	PrintVector(pVector, Size);

	start = omp_get_wtime();
	ParallelResultCalculation(pMatrix, pVector, pResult, Size);
	finish = omp_get_wtime();
	duration = finish - start;


	cout << endl << endl << "Результирующий вектор: ";
	PrintVector(pResult, Size);

	cout << endl << endl << "Запускаем наши тесты для проверки" << endl;
	TestResult(pMatrix, pVector, pResult, Size);
	
	cout << endl << "Затраченое время на выполнения алгоритма: " << duration << endl;

	// Очищаем память и завершаем программу
	ProcessTermination(pMatrix, pVector, pResult);
	return 0;
}

void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int& Size) {
	// Выделяем память
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];

	DataInitializationPreset(pMatrix, pVector, Size);
}

void DataInitializationPreset(double* pMatrix, double* pVector, int Size) {
	pMatrix[0] = 6.0;
	pMatrix[1] = 6.0;
	pMatrix[2] = 5.0;
	pMatrix[3] = 18.0;
	pMatrix[4] = 20.0;

	pMatrix[Size] = 10.0;
	pMatrix[Size + 1] = 9.0;
	pMatrix[Size + 2] = 7.0;
	pMatrix[Size + 3] = 24.0;
	pMatrix[Size + 4] = 30.0;

	pMatrix[2 * Size] = 12.0;
	pMatrix[2 * Size + 1] = 12.0;
	pMatrix[2 * Size + 2] = 13.0;
	pMatrix[2 * Size + 3] = 27.0;
	pMatrix[2 * Size + 4] = 35.0;

	pMatrix[3 * Size] = 8.0;
	pMatrix[3 * Size + 1] = 6.0;
	pMatrix[3 * Size + 2] = 6.0;
	pMatrix[3 * Size + 3] = 15.0;
	pMatrix[3 * Size + 4] = 20.0;

	pMatrix[4 * Size] = 4.0;
	pMatrix[4 * Size + 1] = 5.0;
	pMatrix[4 * Size + 2] = 4.0;
	pMatrix[4 * Size + 3] = 15.0;
	pMatrix[4 * Size + 4] = 15.0;


	pVector[0] = 14.0;
	pVector[1] = 18.0;
	pVector[2] = 32.0;
	pVector[3] = 16.0;
	pVector[4] = 11.0;
}

void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	int i, j;
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			cout << pMatrix[i * RowCount + j] << " ";
		cout << endl;
	}
}

void PrintVector(double* pVector, int Size) {
	for (int i = 0; i < Size; i++)
		cout << pVector[i] << " ";
}

int ParallelFindPivotRow(double* pMatrix, int Size, int Iter) {
	int PivotRow = -1;    // The index of the pivot row
	double MaxValue = 0; // The value of the pivot element
	int i;

	// Выбираем строку где хранится максимальный элемент
#pragma omp parallel
	{
		TThreadPivotRow ThreadPivotRow;
		ThreadPivotRow.MaxValue = 0;
		ThreadPivotRow.PivotRow = -1;
#pragma omp for
		for (i = 0; i < Size; i++) {
			if ((pPivotIter[i] == -1) &&
				(fabs(pMatrix[i * Size + Iter]) > ThreadPivotRow.MaxValue)) {
				ThreadPivotRow.PivotRow = i;
				ThreadPivotRow.MaxValue = fabs(pMatrix[i * Size + Iter]);
			}
		}

#pragma omp critical
		{
			if (ThreadPivotRow.MaxValue > MaxValue) {
				MaxValue = ThreadPivotRow.MaxValue;
				PivotRow = ThreadPivotRow.PivotRow;
			}
		}
	}
	return PivotRow;
}

void TestResult(double* pMatrix, double* pVector, double* pResult, int Size) {
	/* Buffer for storing the vector, that is a result of multiplication
	of the linear system matrix by the vector of unknowns */
	double* pRightPartVector;
	// Flag, that shows wheather the right parts
	// vectors are identical or not
	int equal = 0;
	double Accuracy = 1.e-6; // Comparison accuracy
	pRightPartVector = new double[Size];
	for (int i = 0; i < Size; i++) {
		pRightPartVector[i] = 0;
		for (int j = 0; j < Size; j++) {
			pRightPartVector[i] +=
				pMatrix[i * Size + j] * pResult[j];
		}
	}
	for (int i = 0; i < Size; i++) {
		if (fabs(pRightPartVector[i] - pVector[i]) > Accuracy)
			equal = 1;
	}
	if (equal == 1)
		cout << "Результат параллельного вычисления НЕ правильный!" << endl;
	else
		cout << "Результат параллельного вычисления правильный!" << endl;

	delete[] pRightPartVector;
}

void ParallelBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size) {
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; i--) {
		RowIndex = pPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];

#pragma omp parallel for private (Row)
		for (int j = 0; j < i; j++) {
			Row = pPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}

void ParallelColumnElimination(double* pMatrix, double* pVector, int Pivot, int Iter, int Size) {
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[Pivot * Size + Iter];

#pragma omp parallel for private(PivotFactor) schedule(dynamic,1)
	for (int i = 0; i < Size; i++) {
		if (pPivotIter[i] == -1) {
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++) {
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[Pivot];
		}
	}
}


void ParallelGaussianElimination(double* pMatrix, double* pVector, int Size) {
	int Iter;       // The number of the iteration of the Gaussian
	// elimination
	int PivotRow;   // The number of the current pivot row
	for (Iter = 0; Iter < Size; Iter++) {
		// Finding the pivot row
		PivotRow = ParallelFindPivotRow(pMatrix, Size, Iter); pPivotPos[Iter] = PivotRow;
		pPivotIter[PivotRow] = Iter;
		ParallelColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}

void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {
	// Выделяем память
	pPivotPos = new int[Size];
	pPivotIter = new int[Size];
	for (int i = 0; i < Size; i++) {
		pPivotIter[i] = -1;
	}

	ParallelGaussianElimination(pMatrix, pVector, Size);
	ParallelBackSubstitution(pMatrix, pVector, pResult, Size);

	// Высвобождаем память
	delete[] pPivotPos;
	delete[] pPivotIter;
}

void ProcessTermination(double* pMatrix, double* pVector, double* pResult) {
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
}