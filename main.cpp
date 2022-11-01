#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys, const std::vector<Type> &startPoint,
const std::string &IN_FILE_PATH, Type accuracy = 1e-7){
    // Считываем данные
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);

}

template<typename Type>
void temp_main(){
    std::vector<std::vector<Type>> A; // Матрица левых коэффициентов
    std::vector<Type> b; // Вектор правых коэффициентов
    std::vector<Type> eigList;
    std::vector<std::vector<Type>> Q;
    readMatrix(A, IN_FILE_PATH_3);
    std::cout << findEigenNumsQRMethod(A, eigList, 1e-6, 0) << '\n';
    std::cout << eigList << '\n';
    readMatrix(A, IN_FILE_PATH_3);
    std::cout << findEigenNumsQRMethod(A, eigList, 1e-6, 1) << '\n';
    std::cout << eigList << '\n';
    readMatrix(A, IN_FILE_PATH_3);
    std::cout << findEigenNumsQRMethodHessenberg(A, eigList, 1e-6, 0) << '\n';
    std::cout << eigList << '\n';
    readMatrix(A, IN_FILE_PATH_3);
    std::cout << findEigenNumsQRMethodHessenberg(A, eigList, 1e-6, 1) << '\n';
    std::cout << eigList << '\n';
    
}
int main(){
    temp_main<double>();
    return 0;
}