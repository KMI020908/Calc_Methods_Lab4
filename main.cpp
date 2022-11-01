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
    invertItersMethod(matrix, eigMatrix, eigList, accuracy, 0.1);
    writeEigenVec(eigMatrix, eigList, I_OUT_FILE_PATH);


}

template<typename Type>
void temp_main(){
    std::vector<std::vector<Type>> matrix; 
    std::vector<Type> eigList;
    checkTest(matrix, eigList, IN_FILE_PATH_1, QR_OUT_FILE_PATH_1, I_OUT_FILE_PATH_1);
    checkTest(matrix, eigList, IN_FILE_PATH_2, QR_OUT_FILE_PATH_2, I_OUT_FILE_PATH_2);
    checkTest(matrix, eigList, IN_FILE_PATH_3, QR_OUT_FILE_PATH_3, I_OUT_FILE_PATH_3);
}

int main(){
    temp_main<double>();
    return 0;
}