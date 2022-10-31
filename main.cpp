#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys, const std::vector<Type> &startPoint,
const std::string &IN_FILE_PATH, const std::string &SIMPLE_IT_F_PATH, const std::string &JACOBI_F_PATH, const std::string &RELAXATION_F_PATH, 
Type accuracy = 1e-7, std::vector<double> realSolution = {}){
    // Считываем данные
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);

}

template<typename Type>
void temp_main(){

    std::vector<std::vector<Type>> A; // Матрица левых коэффициентов
    std::vector<Type> b; // Вектор правых коэффициентов
    readData(A, b, IN_FILE_PATH_2);
    std::vector<Type> eigList;
    findEigenNumsQRMethodClassic(A, eigList, 1e-7);
    std::cout << eigList << '\n';


}

int main(){
    temp_main<double>();
    return 0;
}