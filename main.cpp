#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigList,
const std::string &IN_FILE_PATH, const std::string &QR_OUT_FILE_PATH, const std::string &I_OUT_FILE_PATH, Type accuracy = 1e-6, 
bool is3Diag = false){
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
    numOfIters = invertItersMethod(matrix, eigMatrix, eigList, accuracy, false, 0.2);
    writeEigenVec(numOfIters, eigMatrix, eigList, I_OUT_FILE_PATH);
    std::vector<Type> eigVec1;
    std::vector<Type> startVec = { 1.0, 1.0, 1.0 };
    Type lambda1 = invertItersMethodRayleigh(matrix, startVec, eigVec1, accuracy, is3Diag, 1.0);
    std::cout << "Eigen number: " << lambda1 << " ----> " << "Eigen vector: " << eigVec1 << '\n' << '\n';

    //writeEigenVec(numOfIters, startVec, eigVec1, I_OUT_FILE_PATH, true);
}

template<typename Type>
void temp_main(){
    std::vector<std::vector<Type>> matrix; 
    std::vector<Type> eigList;
    Type accuracy = 1e-6;
    //checkTest(matrix, eigList, IN_FILE_PATH_1, QR_OUT_FILE_PATH_1, I_OUT_FILE_PATH_1, accuracy);
    //checkTest(matrix, eigList, IN_FILE_PATH_2, QR_OUT_FILE_PATH_2, I_OUT_FILE_PATH_2, accuracy);
    checkTest(matrix, eigList, IN_FILE_PATH_3, QR_OUT_FILE_PATH_3, I_OUT_FILE_PATH_3, accuracy, true);

    // Бассейн Ньютона
    //readMatrix(matrix, IN_FILE_PATH_3);
    //findEigenNumsQRMethodHessenberg(matrix, eigList);
    //writeNewthonSwPool(matrix, 0.02, I_OUT_FILE_PATH_4, 1e-6, true);

    // Матрица опрератора дифференцирования
    //std::size_t dim = 10;
    //std::size_t dim = 50;
    std::size_t dim = 115;
    Type h = 1e-4;
    Type c = -1 / (h * h);
    std::vector<std::vector<Type>> matrix2;
    matrix2.resize(dim);
    for (size_t i = 0; i < dim; i++){
        matrix2[i].resize(dim, 0.0);
    }
    matrix2[0][0] = -2.0 * c;
    matrix2[0][1] = 1.0 * c;
    matrix2[dim - 1][dim - 1] = -2.0 * c;
    matrix2[dim - 1][dim - 2] = 1.0 * c;
    for (size_t i = 1; i < dim - 1; i++){
        matrix2[i][i] = -2.0 * c;
        matrix2[i][i + 1] = 1.0 * c;
        matrix2[i][i - 1] = matrix2[i][i + 1];
    }
    writeMatrixFile(matrix2, IN_FILE_PATH_5);
    std::vector<Type> startVec(dim, 1.0);
    std::vector<Type> eigVec;
    Type lambda = invertItersMethodRayleigh(matrix2, startVec, eigVec, accuracy, true);
    std::cout << "Eigen number: " << lambda << " ----> " << "Eigen vector: " << eigVec << '\n' << '\n';

    std::vector<std::vector<Type>> eigMatrix;
    std::size_t numOfIters = findEigenNumsQRMethodHessenberg(matrix2, eigList, accuracy, false, true);
    writeEigenData<Type>(numOfIters, eigList, QR_OUT_FILE_PATH_5, true, true, false);
    matrix2[0][0] = -2.0 * c;
    matrix2[0][1] = 1.0 * c;
    matrix2[dim - 1][dim - 1] = -2.0 * c;
    matrix2[dim - 1][dim - 2] = 1.0 * c;
    for (size_t i = 1; i < dim - 1; i++){
        matrix2[i][i] = -2.0 * c;
        matrix2[i][i + 1] = 1.0 * c;
        matrix2[i][i - 1] = matrix2[i][i + 1];
    }
    numOfIters = invertItersMethod(matrix2, eigMatrix, eigList, accuracy, true);
    writeEigenVec(numOfIters, eigMatrix, eigList, I_OUT_FILE_PATH_5);
    writeVectorFile(eigList, IN_FILE_PATH_8);
    writeMatrixFile(eigMatrix, IN_FILE_PATH_9);
}

int main(){
    temp_main<double>();
}