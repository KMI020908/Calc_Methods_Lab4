#ifndef METHODS_H
#define METHODS_H

#include<vector>
#include<cmath>
#include"ioData.h"
#include<iostream>

// Методы
template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7);

template<typename Type>
SOLUTION_FLAG gaussMethodFull(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7); // Метод Гаусса с полным выбором главного элемента

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-7);

template<typename Type>
SOLUTION_FLAG tridiagonalAlgoritm(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, std::vector<Type> &solution);

template<typename Type>
std::size_t simpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> 
&firstVec, std::vector<Type> &solution, Type tao, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-4, std::size_t stopIt = 100000);

template<typename Type>
std::size_t JacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy = 1e-7, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

template<typename Type>
std::size_t relaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy = 1e-7, Type omega = 1.0, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

template<typename Type>
std::size_t relaxationMethodFor3Diag(const std::vector<Type> &a, const std::vector<Type> &b, const std::vector<Type> & c, const std::vector<Type> &d, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy  = 1e-7, Type omega = 1.0, double p = 2.0, Type epsilon_0 = 1e-7, std::size_t stopIt = 100000);

// Вспомогательные функции
template<typename Type>
QUADRATIC_FLAG findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
QUADRATIC_FLAG findQMatrix3Diag(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution); // Найти невязку

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix); // Октаэдрическая норма матрицы

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &matrix); // Число обусловленности с октаэдоической метрикой

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix); // Кубическая норма матрицы

template<typename Type>
Type normOfMatrix(const std::vector<std::vector<Type>> &matrix, double p = 1.0); 

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &matrix); // Число обусловленности с кубической метрикой

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector); // Октаэдрическая норма вектора

template<typename Type>
Type norm2OfVector(const std::vector<Type> &vector); // Квадратичная норма вектора

template<typename Type>
Type normOfVector(const std::vector<Type> &vector, double p = 2.0);

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector); // Кубическая норма вектора

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix,
    SOLUTION_FLAG (*method)(std::vector<std::vector<Type>>& matrix, std::vector<Type>& rVec, std::vector<Type>& sol, Type accuracy) = qrMethod); // Обратная матрица

template<typename Type>
std::size_t transposeMatrix(std::vector<std::vector<Type>> &matrix);

template<typename Type>
Type findLowerBoundOfcond1(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1 = 1e-1, Type delta2 = 1e-3, Type delta3 = 1e-6);

template<typename Type>
Type findLowerBoundOfcondInf(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1 = 1e-1, Type delta2 = 1e-3, Type delta3 = 1e-6);

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(std::vector<std::vector<Type>> &matrix1, const std::vector<std::vector<Type>> &matrix2, std::vector<std::vector<Type>> &result);

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec, std::vector<Type> &result);

// Перегрузки
template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<Type>> &matrix);

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<Type> &vector);

template<typename Type>
std::vector<Type> operator+(const std::vector<Type>& vec1, const std::vector<Type>& vec2);

template<typename Type>
std::vector<Type> operator-(const std::vector<Type>& vec1, const std::vector<Type>& vec2);

template<typename Type>
std::vector<Type> operator*(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec);

template<typename Type>
std::vector<Type> operator*(Type num, const std::vector<Type> &vec);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormSimpleIt(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type tao);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormJacobi(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega = 1.0);

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation2(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega = 1.0);

template<typename Type>
Type findLowerBoundOfIterations(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, const std::vector<Type> &firstVec, Type accuracy, ITERATION_METHOD_FLAG method, 
Type tao, Type omega = 1.0, double p = 2.0);

template<typename Type>
std::size_t findExactItersSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
std::size_t findExactItersJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
std::size_t findExactRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, Type omega, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_SIT(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_JAC(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template<typename Type>
Type findNormOfErrAfterEstIt_REL(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt = SIZE_MAX);

template <typename Type>
std::size_t findEigenNumsQRMethod(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy = 1e-6, bool hasShift = false);

template <typename Type>
QUADRATIC_FLAG getHessenbergMatrix(std::vector<std::vector<Type>> &matrix, Type accuracy = 1e-6, bool isSymmetric = false);

template<typename Type>
std::size_t findEigenNumsQRMethodHessenberg(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy, bool hasShift);

template<typename Type>
std::size_t invertItersMethod(const std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &eigenMatrix, const std::vector<Type> &startEigenList,
Type accuracy = 1e-6, Type omega = 1.0);

template<typename Type>
Type invertItersMethodRayleigh(const std::vector<std::vector<Type>> &matrix, std::vector<Type> &startVec, 
std::vector<Type> &eigenVec, Type accuracy, Type omega);

template<typename Type>
Type dot(const std::vector<Type> &v1, const std::vector<Type> &v2);

#endif