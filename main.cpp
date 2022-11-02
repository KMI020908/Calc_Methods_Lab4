#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigList,
const std::string &IN_FILE_PATH, const std::string &QR_OUT_FILE_PATH, const std::string &I_OUT_FILE_PATH, Type accuracy = 1e-6){
    // QR метод
    std::vector<std::vector<Type>> eigMatrix;
    std::size_t numOfIters = 0;

    readMatrix<Type>(matrix, IN_FILE_PATH);
    numOfIters = findEigenNumsQRMethod<Type>(matrix, eigList, accuracy, 0);
    writeEigenData<Type>(numOfIters, eigList, QR_OUT_FILE_PATH, false, false, false);
    
    readMatrix<Type>(matrix, IN_FILE_PATH);
    numOfIters = findEigenNumsQRMethod<Type>(matrix, eigList, accuracy, 1);
    writeEigenData<Type>(numOfIters, eigList, QR_OUT_FILE_PATH, true, false, true);

    readMatrix<Type>(matrix, IN_FILE_PATH);
    numOfIters = findEigenNumsQRMethodHessenberg<Type>(matrix, eigList, accuracy, 1);
    writeEigenData<Type>(numOfIters, eigList, QR_OUT_FILE_PATH, false, true, true);

    readMatrix<Type>(matrix, IN_FILE_PATH);
    numOfIters = findEigenNumsQRMethodHessenberg<Type>(matrix, eigList, accuracy, 1);
    writeEigenData<Type>(numOfIters, eigList, QR_OUT_FILE_PATH, true, true, true);

    readMatrix<Type>(matrix, IN_FILE_PATH);
    numOfIters = invertItersMethod(matrix, eigMatrix, eigList, accuracy, 0.2);
    writeEigenVec(numOfIters, eigMatrix, eigList, I_OUT_FILE_PATH);
    std::vector<Type> eigVec1;
    std::vector<Type> startVec = { 0.777976, 0.582499, -0.229449, 0.0529363 };
    Type lambda1 = invertItersMethodRayleigh(matrix, startVec, eigVec1, accuracy, 0.2);
    std::cout << "Eigen number: " << lambda1 << " ----> " << "Eigen vector: " << eigVec1 << '\n' << '\n';

    //writeEigenVec(numOfIters, startVec, eigVec1, I_OUT_FILE_PATH, true);
}

template<typename Type>
void temp_main(){
    std::vector<std::vector<Type>> matrix; 
    std::vector<Type> eigList;
    checkTest(matrix, eigList, IN_FILE_PATH_1, QR_OUT_FILE_PATH_1, I_OUT_FILE_PATH_1);
    //checkTest(matrix, eigList, IN_FILE_PATH_2, QR_OUT_FILE_PATH_2, I_OUT_FILE_PATH_2);
    //checkTest(matrix, eigList, IN_FILE_PATH_3, QR_OUT_FILE_PATH_3, I_OUT_FILE_PATH_3);
}

int main(){
    temp_main<double>();
    /*
    std::vector<std::vector<double>> A;
    readMatrix(A, IN_FILE_PATH_3);
    getHessenbergMatrix(A);
    std::cout << A;
    */
    return 0;
}